import os
from ..tools import decorators
from ..tools import classes
from ..tools.logger import logger
from ..tools.script_tools import create_folder


@decorators.user_choice
@decorators.track_run
def rdkit_sampling(
                        receptor_obj,
                        ligase_obj,
                        protac_obj,
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
    from rdkit.Geometry.rdGeometry import Point3D
    from func_timeout import func_timeout
    from ..tools.structure_tools import smiles2smarts
    from ..run.megadock import rotate_atoms

    # make folder where the linkers for all pose objs will be stored
    create_folder(protac_poses_folder)

    # open receptor ligand - will not change
    reclig = Chem.MolFromMol2File(receptor_obj.lig_file)
    # get initial rotation information from ligase obj
    ref_rotation  = ligase_obj.rotate


    pose_objs = ligase_obj.active_confs()
    for pose_obj in pose_objs:

        # make folder for each pose obj linker file to be saved
        linker_folder = os.path.join(protac_poses_folder, f'protein_pose_{pose_obj.pose_number}')
        create_folder(linker_folder)

        # open ligase ligand - changes every loop
        liglig = Chem.MolFromMol2File(ligase_obj.lig_file)
        conf = liglig.GetConformer()

        # get final rotation information for the pose
        pose_rotation = pose_obj.rotate
        
        # rotate ligase ligand to new position
        for i in range(liglig.GetNumAtoms()):
            x, y, z = conf.GetAtomPosition(i)
            newX, newY, newZ = rotate_atoms((x, y, z), ref_rotation=ref_rotation, pose_rotation=pose_rotation)
            conf.SetAtomPosition(i,Point3D(newX, newY, newZ))

        # combine both ligs
        reference_ligs = Chem.CombineMols(liglig, reclig)

        #  convert and get smiles
        pose_lig_smiles = Chem.MolToSmiles(liglig)
        receptor_lig_smiles = Chem.MolToSmiles(reclig)
        
        # for each of their smiles, transform into wildcard smarts
        receptor_lig_smarts_ = smiles2smarts(receptor_lig_smiles)
        pose_lig_smarts_ = smiles2smarts(pose_lig_smiles)
        receptor_lig_smarts = Chem.MolFromSmarts(receptor_lig_smarts_)
        pose_lig_smarts = Chem.MolFromSmarts(pose_lig_smarts_)
        
        # open protac smiles and addhs
        protac = Chem.MolFromSmiles(protac_obj.smiles)
        protac = Chem.AddHs(protac)

        # find atoms in the protac that match the wildcards:
        receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]
        # find atoms in the reference structure that match the wildcards,
        # they are going to give the coordinates we have to match:
        receptor_lig_coords = reference_ligs.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_coords = reference_ligs.GetSubstructMatches(pose_lig_smarts)[0]

        # make a new protac object and save its ligand atom indices if not done before
        if protac_obj.index_ligs == None:
            indices = list(receptor_lig_indices)
            indices.extend(list(pose_lig_indices))
            protac_obj.index_ligs = indices

        # build the dict that will contain the mapping btwn indices and coords:
        coordmap = {}
        mol = reference_ligs.GetConformer()
        for i in range(len(receptor_lig_indices)):
            atom_ix = receptor_lig_indices[i]
            atom_coord = receptor_lig_coords[i]
            x, y, z = mol.GetAtomPosition(atom_coord)
            coordmap[atom_ix] = Point3D(x, y, z)
        for i in range(len(pose_lig_indices)):
            atom_ix = pose_lig_indices[i]
            atom_coord = pose_lig_coords[i]
            x, y, z = mol.GetAtomPosition(atom_coord)
            coordmap[atom_ix] = Point3D(x, y, z)

        # make ProtacPose obj:
        protac_pose_obj = classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)
        # sample its conformations!
        kwargs = {
            'mol':protac, 'coordMap':coordmap,
            'numConfs':rdkit_number_of_confs,
            'enforceChirality':False
        }
        try:
            func_timeout(time_tolerance, AllChem.EmbedMultipleConfs, kwargs=kwargs)
        except:
            logger.warning(f"rdkit timeout for pose {pose_obj.pose_number}")
            pass

        # make pairlist to align the protac confs to the ref struct
        pairs = [(receptor_lig_indices[i], receptor_lig_coords[i]) for i in range(len(receptor_lig_indices))]
        pairs.extend([(pose_lig_indices[i], pose_lig_coords[i]) for i in range(len(pose_lig_indices))])
        
        # for each conformation, align and write to single sdf file,
        # capturing only the poses that obey the rmsd tolerance
        protac_file = os.path.join(linker_folder, 'protac_embedded_confs.sdf')
        try:
            with open(protac_file, 'a+') as confs_file:
                for i in range(rdkit_number_of_confs):
                    # make linker conf object
                    linker_conf = classes.LinkerConf(parent=protac_pose_obj, conf_number=i)
                    rmsd = Chem.rdMolAlign.AlignMol(protac, reference_ligs, atomMap=pairs, prbCid=i)
                    if rmsd <= rmsd_tolerance:
                        molblock = Chem.MolToMolBlock(protac, confId=i, kekulize=False)
                        confs_file.write(f'conf {i}')
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

        # update pose_obj with new linker file
        pose_obj.protac_file = protac_file


@decorators.track_run
@decorators.user_choice
def dock6_score(pose_objs, dock6_root):
    """
    Use dock6 continuous score to rank the linker conformations.
    If no linker conf has negative energy, then the protein pose is not viable
    """
    import subprocess
    from ..definitions import ROOT_DIR

    logger.info('Scoring protacs coformations using dock6 continuous score')

    for pose_obj in pose_objs:

        # use obabel to convert sdf to mol2:
        mol2_file = pose_obj.protac_file.replace('sdf','mol2')
        command = [
            'obabel', '-isdf', pose_obj.protac_file,
            '-omol2', '-O', mol2_file
        ]
        subprocess.run(command)
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
                found = [int(i.split(' ')[-1]) for i in found]
            else:
                found = [float(i) for i in found]
            data[key] = found
        file_text = None

        return(data)
 
    for pose_obj in pose_objs:

        folder_path = os.path.dirname(pose_obj.protac_file)
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


# def detect_clashes(receptor_obj, pose_objs, clash_tolerance):
#     """
#     Use biopython to detect clashes between proteins and linker poses
#     """
#     # pose_objs where linker_gen is True

#     receptor_struct = receptor_obj.get_protein_struct()

#     for pose_obj in pose_objs:
