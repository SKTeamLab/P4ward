from ..tools.logger import logger
from ..tools import decorators
from pathlib import Path
import subprocess


def load_biopython_structures(structure_file, mol2=False):
    """
    load a structure as biopython object. If file in mol2 format,
    convert first using obabel. Intended for use with the ligands.
    """

    from Bio.PDB.PDBParser import PDBParser as pdbp
    parser = pdbp(PERMISSIVE=1)

    if mol2:
        from io import StringIO

        run = subprocess.run(
            ['obabel','-imol2',structure_file,'-opdb'],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True
        )
        pdb_block = run.stdout
        pdbfile = StringIO(pdb_block)

        structure_obj = parser.get_structure('structure', pdbfile)

    else:
        structure_obj = parser.get_structure('structure', structure_file)

    return(structure_obj)


@decorators.user_choice
@decorators.track_run
def get_protac_dist_cuttoff(
        protac_objs,
        reclig_file,
        liglig_file,
        dist_cutoff,
        sampling_type,
        unbound_protac_num_confs
):

    if dist_cutoff == 'auto':

        logger.info("Ligands distance cutoff set to automatic.")

        for protac_obj in protac_objs:
            logger.info(f'Sampling unbound conformations for {protac_obj.name} to determine distance cutoff.')

            import numpy as np
            from rdkit import Chem
            from rdkit.Chem import rdFMCS, AllChem
            
            def center_of_mass(atom_indices, conf):

                total_mass = 0.0
                xs = []; ys = []; zs = []
                for j in atom_indices:
                    atom = conf.GetAtomWithIdx(j)
                    mass = atom.GetMass()
                    total_mass += mass
                    x,y,z = conf.GetConformer().GetAtomPosition(j)
                    xs.append(mass*x)
                    ys.append(mass*y)
                    zs.append(mass*z)
                center = np.array([np.sum(xs)/total_mass, np.sum(ys)/total_mass, np.sum(zs)/total_mass])
                return(center)
            
            if sampling_type == '2D':

                unbound_protac_num_confs = 1
                protac = Chem.MolFromSmiles(protac_obj.smiles)
                AllChem.Compute2DCoords(protac)
                num_sampled_confs = 1
            
            elif sampling_type == '3D':

                protac_obj.sample_unbound_confs(num_unbound_confs=unbound_protac_num_confs)
                protac = protac_obj.unbound_confs
                num_sampled_confs = protac_obj.num_confs

                if num_sampled_confs == 0 or num_sampled_confs <= 0.1 * unbound_protac_num_confs:
                    logger.warning(
                        f"RDKit sampled {num_sampled_confs} but {unbound_protac_num_confs} were requested. "+
                        f"This could be a challenging molecule, consider increasing unbound_protac_num_confs and rdkit_number_of_confs."
                    )

            reclig = Chem.MolFromMol2File(str(reclig_file), sanitize=False, cleanupSubstructures=False)
            liglig = Chem.MolFromMol2File(str(liglig_file), sanitize=False, cleanupSubstructures=False)

            receptor_lig_smarts_ = rdFMCS.FindMCS([protac, reclig]).smartsString
            pose_lig_smarts_ = rdFMCS.FindMCS([protac, liglig]).smartsString
            receptor_lig_smarts = Chem.MolFromSmarts(receptor_lig_smarts_)
            pose_lig_smarts = Chem.MolFromSmarts(pose_lig_smarts_)
            receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
            pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]
        
            distances = []

            # iterate through how many conformations were actually sampled for the protac
            for i in range(num_sampled_confs):

                conf = Chem.Mol(protac, confId=i)
                center_reclig = center_of_mass(receptor_lig_indices, conf)
                center_liglig = center_of_mass(pose_lig_indices, conf)

                distance = np.linalg.norm(center_liglig - center_reclig)
                distances.append(distance)

            cutoff = np.mean(distances)
            protac_obj.dist_cutoff = cutoff

            logger.info(f"Setting distance cutoff to {cutoff}")
    
    else:
        for protac_obj in protac_objs:
            protac_obj.dist_cutoff = float(dist_cutoff)
        logger.info(f"Ligands distance cutoff set to {dist_cutoff}.")


def structure_proximity(struct1, struct2, dist_cutoff=None):
    """
    Calculate the distance between the center of mass of two biopython struct objects
    and check if it is <= than a distance cutoff.
    """

    import numpy as np

    distance = np.linalg.norm(
        struct1.center_of_mass() - struct2.center_of_mass()
    )
    return(distance, distance <= dist_cutoff)


def get_rmsd(obj1, obj2, fit=False, ca=False, backbone=True):
    """
    Use biopython to get the atomic coordinates of both objects
    and the pip tool rmsd.rmsd to calulate the rmsd. Returns rmsd value
    """

    import numpy as np
    import rmsd

    objects = {
        'obj1': obj1,
        'obj2': obj2,
        'obj1_coords':[],
        'obj2_coords':[]
    }
    backbone_atoms = ['C','CA','N','O']
    
    for obj in ('obj1', 'obj2'):
        for res in objects[obj].get_residues():
            for atom in res.get_atoms():
                if ca:
                    if atom.get_name() == 'CA':
                        objects[f'{obj}_coords'].append(atom.get_coord())
                    else: continue
                elif backbone:
                    if atom.get_name() in backbone_atoms:
                        objects[f'{obj}_coords'].append(atom.get_coord())
                    else: continue
                else:
                    objects[f'{obj}_coords'].append(atom.get_coord())
    
    obj1_coords = np.asarray(objects['obj1_coords'])
    obj2_coords = np.asarray(objects['obj2_coords'])

    if fit:
        rmsd = rmsd.kabsch_rmsd(obj1_coords, obj2_coords, translate=True)
    else:
        rmsd = rmsd.rmsd(obj1_coords, obj2_coords)
    return(rmsd)


def pymol_combine(*args, out_filename='combined.pdb', assign_chains=True):

    import string
    import pymol2

    with pymol2.PyMOL() as pm:

        basenames = [Path(filename).stem for filename in args]
        for filename in args:
            pm.cmd.load(filename)
        
        # loop through models and assign a separate chain to each
        if assign_chains:
            chain_letters = list(string.ascii_lowercase)
            for i in range(len(basenames)):
                model = basenames[i]
                chain = chain_letters[i]
                pm.cmd.alter(model, f'chain="{chain}"')

        ## combine all and save
        pm.cmd.create('combined', ' '.join(basenames))
        pm.cmd.save(out_filename, 'combined')


def pymol_align(target_file, moving_file, outfilename):

    import pymol2

    with pymol2.PyMOL() as pm:

        pm.cmd.load(target_file, 'target')
        pm.cmd.load(moving_file, 'moving')
        pm.cmd.align('moving', 'target')
        pm.cmd.save(outfilename, 'moving')

    return(outfilename)


def write_charges(*charge_lists, filepath):

    charges = []
    for charge_list in charge_lists:
        charges.extend(charge_list)
    
    rawfile = open(filepath, 'r').read().splitlines()

    for coords, charge in charges:
        
        coord_str = f"{format(coords[0],'.3f')}\t{format(coords[1],'.3f')}\t{format(coords[2],'.3f')}\t"
        for i in range(len(rawfile)):
            line = rawfile[i]
            if coord_str in line:
                line = line.replace('0.000', str(charge)) # format(i,'.3f')
                rawfile[i] = line
    
    with open(filepath, 'w') as newfile:
        newfile.write('\n'.join(rawfile))


def get_coords_array(pose_objs, ligase_obj):

    import numpy as np
    from ..run.megadock import rotate_atoms

    a,c,d = ligase_obj.get_triad_points()

    coords = []

    for pose_obj in pose_objs:

        a_rot = rotate_atoms(tuple(a), ref_rotation=ligase_obj.rotate, pose_rotation=pose_obj.rotate)
        c_rot = rotate_atoms(tuple(c), ref_rotation=ligase_obj.rotate, pose_rotation=pose_obj.rotate)
        d_rot = rotate_atoms(tuple(d), ref_rotation=ligase_obj.rotate, pose_rotation=pose_obj.rotate)

        coords.append([*a_rot, *c_rot, *d_rot])

    coords = np.asarray(coords)
    return(coords)
