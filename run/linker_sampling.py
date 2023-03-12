import os
from ..tools import decorators
from ..tools.logger import logger
from ..tools.script_tools import create_folder


@decorators.user_choice
@decorators.track_run
def rdkit_sampling(
                        receptor_obj,
                        ligase_obj,
                        protac,
                        rdkit_number_of_confs,
                        protac_poses_folder,
                        rmsd_tolerance
):

    """
    use rdkit to sample linker conformations while keeping ligands restrained
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry.rdGeometry import Point3D
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
        protac_ = open('protac.smiles', 'r').read()
        protac = Chem.MolFromSmiles(protac_)
        protac = Chem.AddHs(protac)

        # find atoms in the protac that match the wildcards:
        receptor_lig_indices = protac.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_indices = protac.GetSubstructMatches(pose_lig_smarts)[0]
        # find atoms in the reference structure that match the wildcards,
        # they are going to give the coordinates we have to match:
        receptor_lig_coords = reference_ligs.GetSubstructMatches(receptor_lig_smarts)[0]
        pose_lig_coords = reference_ligs.GetSubstructMatches(pose_lig_smarts)[0]

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

        # sample the conformations!
        AllChem.EmbedMultipleConfs(protac, coordMap=coordmap, numConfs=rdkit_number_of_confs)

        # make pairlist to align the protac confs to the ref struct
        pairs = [(receptor_lig_indices[i], receptor_lig_coords[i]) for i in range(len(receptor_lig_indices))]
        pairs.extend([(pose_lig_indices[i], pose_lig_coords[i]) for i in range(len(pose_lig_indices))])
        
        # for each conformation, align and write to single sdf file,
        # capturing only the poses that obey the rmsd tolerance
        protac_file = os.path.join(linker_folder, 'protac_embedded_confs.sdf')
        try:
            with open(protac_file, 'a+') as confs_file:
                for i in range(rdkit_number_of_confs):
                    rmsd = Chem.rdMolAlign.AlignMol(protac, reference_ligs, atomMap=pairs, prbCid=i)
                    if rmsd <= rmsd_tolerance:
                        molblock = Chem.MolToMolBlock(protac, confId=i, kekulize=False)
                        confs_file.write(f'conf {i}')
                        confs_file.write(molblock)
                        confs_file.write('$$$$\n')
        except:
            logger.info(f'No conformation possible for pose {pose_obj.pose_number}')

        # update pose_obj with new linker file
        pose_obj.protac_file = protac_file
