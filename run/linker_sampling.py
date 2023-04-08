import os
import copy
import numpy as np
from ..tools import decorators
from ..tools import classes
from ..tools.logger import logger
from ..tools.script_tools import create_folder


@decorators.user_choice
#@decorators.track_run
def rdkit_sampling(
                        receptor_obj,
                        ligase_obj,
                        protac_obj,
                        extend_flexible_small_linker,
                        neighbour_number,
                        min_linker_length,
                        rdkit_number_of_confs,
                        protac_poses_folder,
                        rmsd_tolerance,
                        time_tolerance
):

    """
    use rdkit to sample linker conformations while keeping ligands restrained
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import rdFMCS
    from rdkit.Geometry.rdGeometry import Point3D
    from func_timeout import func_timeout
    from ..run.megadock import rotate_atoms


    def make_indices(receptor_lig_indices, pose_lig_indices):

        indices_ligs = list(receptor_lig_indices)
        indices_ligs.extend(list(pose_lig_indices))

        indices_link = []
        for atom in Chem.RemoveHs(protac).GetAtoms():
            ix = atom.GetIdx()
            if ix not in indices_ligs:
                indices_link.append(ix)
        
        return(indices_ligs, indices_link)
    

    def get_matches(protac, reclig, liglig, reference_ligs):

        receptor_lig_smarts_ = rdFMCS.FindMCS([protac, reclig]).smartsString
        pose_lig_smarts_ = rdFMCS.FindMCS([protac, liglig]).smartsString
        receptor_lig_smarts = Chem.MolFromSmarts(receptor_lig_smarts_)
        pose_lig_smarts = Chem.MolFromSmarts(pose_lig_smarts_)
        
        # find atoms in the protac that match the wildcards:
        receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]
        # find atoms in the reference structure that match the wildcards,
        # they are going to give the coordinates we have to match:
        receptor_lig_coords = reference_ligs.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_coords = reference_ligs.GetSubstructMatches(pose_lig_smarts)[0]

        matches = {
            'receptor_lig_indices' : list(receptor_lig_indices),
            'pose_lig_indices'     : list(pose_lig_indices),
            'receptor_lig_coords'  : list(receptor_lig_coords),
            'pose_lig_coords'      : list(pose_lig_coords)
        }
        return(matches)


    def check_linker_size(protac, indices_link, min_linker_length, neighbour_number, matches):

        if len(indices_link) <= min_linker_length:
            logger.warning(f"Linker length is {len(indices_link)}, extending flexibility to neighbouring atoms.")

            for _ in range(neighbour_number):
            
                neighbours = []
                protac_noh = Chem.RemoveHs(protac)
                for linker_atom in indices_link:                    # for each linker atom:
                    atom = protac_noh.GetAtomWithIdx(linker_atom)   # get the atom obj
                    nbs = [i.GetIdx() for i in atom.GetNeighbors()] # get idx of the atom's neighbours (nbs)
                    nbs = [i for i in nbs if i not in indices_link] # only keep nbs that are not in the linker already
                    neighbours.extend(nbs)                          # add them to neighbour list
                neighbours = list(np.unique(neighbours))            # remove duplicates

                receptor_lig_indices = dict(zip(range(len(matches['receptor_lig_indices'])), matches['receptor_lig_indices']))
                receptor_lig_coords = dict(zip(range(len(matches['receptor_lig_coords'])), matches['receptor_lig_coords']))            
                for i in range(len(receptor_lig_indices)):    # remove nbs from receptor_lig_indices and coords
                    if receptor_lig_indices[i] in neighbours:
                        receptor_lig_indices.pop(i)
                        receptor_lig_coords.pop(i)
                matches['receptor_lig_indices'] = list(receptor_lig_indices.values())
                matches['receptor_lig_coords'] = list(receptor_lig_coords.values())

                pose_lig_indices = dict(zip(range(len(matches['pose_lig_indices'])), matches['pose_lig_indices']))
                pose_lig_coords = dict(zip(range(len(matches['pose_lig_coords'])), matches['pose_lig_coords']))
                for i in range(len(pose_lig_indices)):        # remove nbs from pose_lig_indices and coords
                    if pose_lig_indices[i] in neighbours:
                        pose_lig_indices.pop(i)
                        pose_lig_coords.pop(i)
                matches['pose_lig_indices'] = list(pose_lig_indices.values())
                matches['pose_lig_coords'] = list(pose_lig_coords.values())

                _, indices_link = make_indices(matches['receptor_lig_indices'], matches['pose_lig_indices'])
            
            logger.info(f"There are now {len(indices_link)} atoms treated as flexible.")
        
        else:
            logger.info(f"Linker length is {len(indices_link)}, not extending flexibility.")
        
        return(matches)
           

    # make folder where the linkers for all pose objs will be stored
    create_folder(protac_poses_folder)
    # open receptor ligand - will not change, get its rotation info
    reclig = Chem.MolFromMol2File(receptor_obj.lig_file, sanitize=False, cleanupSubstructures=False)
    ref_rotation  = ligase_obj.rotate
    # open ligase ligand - raw initial position
    liglig = Chem.MolFromMol2File(ligase_obj.lig_file, sanitize=False, cleanupSubstructures=False)
    # combine reclig and liglig into a single molecule
    reference_ligs = Chem.CombineMols(liglig, reclig)
    # open protac and addHs
    protac = Chem.MolFromSmiles(protac_obj.smiles)
    protac = Chem.AddHs(protac)

    matches = get_matches(protac, reclig, liglig, reference_ligs)
    indices_ligs, indices_link = make_indices(matches['receptor_lig_indices'], matches['pose_lig_indices'])
    alignment_pairs = [(matches['receptor_lig_indices'][i], matches['receptor_lig_coords'][i]) for i in range(len(matches['receptor_lig_indices']))]
    alignment_pairs.extend([(matches['pose_lig_indices'][i], matches['pose_lig_coords'][i]) for i in range(len(matches['pose_lig_indices']))])


    if extend_flexible_small_linker:
        matches = check_linker_size(protac, indices_link, min_linker_length, neighbour_number, matches)
        indices_ligs, indices_link = make_indices(matches['receptor_lig_indices'], matches['pose_lig_indices'])
        alignment_pairs = [(matches['receptor_lig_indices'][i], matches['receptor_lig_coords'][i]) for i in range(len(matches['receptor_lig_indices']))]
        alignment_pairs.extend([(matches['pose_lig_indices'][i], matches['pose_lig_coords'][i]) for i in range(len(matches['pose_lig_indices']))])

    # save protac_obj ligand and linker atom indices if not done before
    if protac_obj.indices_ligs == None:
        protac_obj.indices_ligs = indices_ligs
        protac_obj.indices_link = indices_link


    pose_objs = ligase_obj.active_confs()
    for pose_obj in pose_objs:

        # make folder for each pose obj linker file to be saved
        linker_folder = os.path.join(protac_poses_folder, f'protein_pose_{pose_obj.pose_number}')
        create_folder(linker_folder)

        # open ligase ligand - changes every loop
        liglig_rotate = copy.deepcopy(liglig)
        conf = liglig_rotate.GetConformer()

        # get final rotation information for the pose
        pose_rotation = pose_obj.rotate
        
        # rotate ligase ligand to new position
        for i in range(liglig_rotate.GetNumAtoms()):
            x, y, z = conf.GetAtomPosition(i)
            newX, newY, newZ = rotate_atoms((x, y, z), ref_rotation=ref_rotation, pose_rotation=pose_rotation)
            conf.SetAtomPosition(i,Point3D(newX, newY, newZ))

        # combine both ligs
        reference_ligs_rotate = Chem.CombineMols(liglig_rotate, reclig)

        # copy protac to be embedded
        protac_embed = copy.deepcopy(protac)

        # build the dict that will contain the mapping btwn indices and coords:
        coordmap = {}
        mol = reference_ligs_rotate.GetConformer()
        for i in range(len(matches['receptor_lig_indices'])):
            atom_ix = matches['receptor_lig_indices'][i]
            atom_coord = matches['receptor_lig_coords'][i]
            x, y, z = mol.GetAtomPosition(atom_coord)
            coordmap[atom_ix] = Point3D(x, y, z)
        for i in range(len(matches['pose_lig_indices'])):
            atom_ix = matches['pose_lig_indices'][i]
            atom_coord = matches['pose_lig_coords'][i]
            x, y, z = mol.GetAtomPosition(atom_coord)
            coordmap[atom_ix] = Point3D(x, y, z)

        # make ProtacPose obj:
        protac_pose_obj = classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)
        # sample its conformations!
        kwargs = {
            'mol':protac_embed, 'coordMap':coordmap,
            'numConfs':rdkit_number_of_confs,
            'enforceChirality':False
        }
        try:
            func_timeout(time_tolerance, AllChem.EmbedMultipleConfs, kwargs=kwargs)
        except:
            logger.warning(f"rdkit timeout for pose {pose_obj.pose_number}")
            pass

        # for each conformation, align and write to single sdf file,
        # capturing only the poses that obey the rmsd tolerance
        protac_file = os.path.join(linker_folder, 'protac_embedded_confs.sdf')
        try:
            with open(protac_file, 'a+') as confs_file:
                for i in range(rdkit_number_of_confs):
                    # make linker conf object
                    linker_conf = classes.LinkerConf(parent=protac_pose_obj, conf_number=i)
                    rmsd = Chem.rdMolAlign.AlignMol(protac_embed, reference_ligs_rotate, atomMap=alignment_pairs, prbCid=i)
                    logger.debug(f"pose: {pose_obj.pose_number}, conf: {i}, rmsd: {rmsd}")
                    if rmsd <= rmsd_tolerance:
                        molblock = Chem.MolToMolBlock(protac_embed, confId=i, kekulize=False)
                        confs_file.write(f'conf_{i}')
                        confs_file.write(molblock)
                        confs_file.write('$$$$\n')
                        linker_conf.active = True
                    else:
                        linker_conf.active = False
            protac_pose_obj.active = True
            protac_pose_obj.file = protac_file
        except:
            logger.info(f'No conformation possible for pose {pose_obj.pose_number}')
            protac_pose_obj.active = False



@decorators.track_run
@decorators.user_choice
def dock6_score(pose_objs, dock6_root, linkers_only):
    """
    Use dock6 continuous score to rank the linker conformations.
    If no linker conf has negative energy, then the protein pose is not viable
    """
    import subprocess
    from ..definitions import ROOT_DIR
    from ..tools.structure_tools import obabel_convert

    logger.info('Scoring protacs coformations using dock6 continuous score')

    for pose_obj in pose_objs:

        # use obabel to convert sdf to mol2:
        mol2_file = obabel_convert(pose_obj.protac_pose.file, 'sdf', 'mol2')
        logger.info('Converted conformations sdf file into mol2')
    
        folder_path = os.path.dirname(mol2_file)

        # use chimerax to generate complex mol2 file:
        complex_mol2_file = os.path.join(folder_path, 'complex.mol2')
        command = (
             f"open {pose_obj.parent.file};"
            +f"open {pose_obj.file};"
            +f"combine modelId 9;"
            +f"del #1,2;"
            +f"dockprep acMethod gasteiger;"
            +f"save {complex_mol2_file} models #9"
        )
        subprocess.run(['chimerax', '--nogui'], input=command, encoding='ascii')
        logger.info('Generated combined protein complex mol2 file.')

        # prepare input file
        dock6_input = open(os.path.join(ROOT_DIR, 'inputs', 'dock_score_protacs.in')).read()
        dock6_input = dock6_input.replace('[complex_prep.mol2]',complex_mol2_file)
        dock6_input = dock6_input.replace('[protac_file.mol2]', mol2_file)
        dock6_input = dock6_input.replace('[out_folder]', folder_path)
        dock6_input = dock6_input.replace('[dock6_root]', dock6_root)

        with open(os.path.join(folder_path, 'dock.in'), 'w+') as dock_file:
            dock_file.write(dock6_input)
        
        command = [
            os.path.join(dock6_root, 'bin', 'dock6'),
            '-i', os.path.join(folder_path, 'dock.in'),
            '-o', os.path.join(folder_path, 'dock.out')
        ]
        subprocess.run(command)
        logger.info(f'Ran dock6 rescore for protac conformations on pose {pose_obj.pose_number}')


@decorators.track_run
@decorators.user_choice
def capture_dock6_scores(pose_objs, filter_linkers):

    def parse_file(file_path):
        import re

        data = {}
        patterns = {
            'conf_number':'Name:\s*(.*)\n',
            'score':'Continuous_Score:\s*(.*)\n',
            'vdw':'Continuous_vdw_energy:\s*(.*)\n',
            'es':'Continuous_es_energy:\s*(.*)\n'
        }
        file_text = open(file_path,'r').read()
        
        for key in patterns:
            found = re.findall(patterns[key], file_text)
            if key == 'conf_number':
                found = [int(i.split('_')[-1]) for i in found]
            else:
                found = [float(i) for i in found]
            data[key] = found
        file_text = None

        return(data)
 
    for pose_obj in pose_objs:

        folder_path = os.path.dirname(pose_obj.protac_pose.file)
        scored_file = 'protac_scored.mol2'

        data = parse_file(os.path.join(folder_path, scored_file))
        pose_obj.linker_scores = data

        if filter_linkers:
            active_linkers = [i<0 for i in data['score']]
            if any(active_linkers):
                active_linkers_ix = []
                for i in range(len(active_linkers)):
                    if active_linkers[i]:
                        active_linkers_ix.append(data['conf_number'][i])
                pose_obj.active_linkers = active_linkers_ix
            else:
                pose_obj.active_linkers = None
                pose_obj.active = False

@decorators.track_run
@decorators.user_choice
def detect_clashes(
        receptor_obj,
        protac_obj,
        pose_objs,
        protac_poses_folder,
        clash_threshold,
        restrict_clash_to_linker,
        filter_clashed,
        max_clashes_allowed
):
    """
    Use biopython to detect clashes between proteins and linker poses
    """

    from Bio.PDB import Selection, NeighborSearch
    from ..tools.structure_tools import obabel_convert, load_biopython_structures
    from ..tools.script_tools import create_folder
    
    receptor_struct = receptor_obj.get_protein_struct()

    # only looop the pose_objs if their protac_pose is active
    pose_objs = [i for i in pose_objs if i.protac_pose.active]
    for pose_obj in pose_objs:
        pose_struct = pose_obj.get_rotated_struct(struct_type='protein')

        # convert their poses and split
        split_confs_folder = os.path.join(protac_poses_folder, f'protein_pose_{pose_obj.pose_number}', 'split_confs')
        create_folder(split_confs_folder)
        obabel_convert(pose_obj.protac_pose.file, 'sdf', 'pdb', split=True, split_folder=split_confs_folder)

        # initialize neighbour search
        proteins = Selection.unfold_entities(receptor_struct, 'A')
        proteins.extend(Selection.unfold_entities(pose_struct, 'A'))
        ns = NeighborSearch(proteins)

        # only loop the linker_confs that are active
        linker_confs = [i for i in pose_obj.protac_pose.linker_confs if i.active]
        for linker_conf in linker_confs:

            conf_file = f'conf_{linker_conf.conf_number}.pdb'
            conf_struct = load_biopython_structures(os.path.join(split_confs_folder, conf_file))

            atom_selection = Selection.unfold_entities(conf_struct, 'A')
            if restrict_clash_to_linker:
                atom_selection = [i for i in atom_selection if i not in protac_obj.indices_ligs]

            clash_count = 0
            for atom in atom_selection:
                close_atoms = ns.search(atom.coord, clash_threshold)
                if len(close_atoms) > 0:
                    clash_count+=1
            
            linker_conf.clash_count = clash_count
            if filter_clashed and clash_count > max_clashes_allowed:
                linker_conf.active = False



