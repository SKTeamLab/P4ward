import copy
from rdkit.Geometry.rdGeometry import Point3D
from rdkit import Chem
from rdkit.Chem import AllChem
from func_timeout import func_timeout
# from ..tools.logger import logger
from ..tools import classes
from ..run.megadock import rotate_atoms


def conf_sampling(params, pose_obj, protac_obj, logger):

    # make folder for each pose obj linker file to be saved
    linker_folder = params['protac_poses_folder'] /f'protac_{protac_obj.name}' / f'protein_pose_{pose_obj.pose_number}'
    linker_folder.mkdir(exist_ok=True, parents=True)

    # open ligase ligand - changes every loop
    liglig_rotate = copy.deepcopy(protac_obj.parameters['liglig'])
    conf = liglig_rotate.GetConformer()

    # get final rotation information for the pose
    pose_rotation = pose_obj.rotate

    # rotate ligase ligand to new position
    for i in range(liglig_rotate.GetNumAtoms()):
        x, y, z = conf.GetAtomPosition(i)
        newX, newY, newZ = rotate_atoms((x, y, z), ref_rotation=protac_obj.parameters['ref_rotation'], pose_rotation=pose_rotation)
        conf.SetAtomPosition(i,Point3D(newX, newY, newZ))

    # combine both ligs
    reference_ligs_rotate = Chem.CombineMols(liglig_rotate, protac_obj.parameters['reclig'])

    # copy protac to be embedded
    protac_embed = copy.deepcopy(protac_obj.parameters['protac'])

    # build the dict that will contain the mapping btwn indices and coords:
    coordmap = {}
    mol = reference_ligs_rotate.GetConformer()
    for i in range(len(protac_obj.parameters['matches']['receptor_lig_indices'])):
        atom_ix = protac_obj.parameters['matches']['receptor_lig_indices'][i]
        atom_coord = protac_obj.parameters['matches']['receptor_lig_coords'][i]
        x, y, z = mol.GetAtomPosition(atom_coord)
        coordmap[atom_ix] = Point3D(x, y, z)
    for i in range(len(protac_obj.parameters['matches']['pose_lig_indices'])):
        atom_ix = protac_obj.parameters['matches']['pose_lig_indices'][i]
        atom_coord = protac_obj.parameters['matches']['pose_lig_coords'][i]
        x, y, z = mol.GetAtomPosition(atom_coord)
        coordmap[atom_ix] = Point3D(x, y, z)
    
    # sample protac conformations!
    kwargs = {
        'mol':protac_embed, 'coordMap':coordmap,
        'numConfs':params['rdkit_number_of_confs'],
        'randomSeed':params['rdkit_random_seed'],
        'enforceChirality':False
    }
    try:
        func_timeout(params['time_tolerance'], AllChem.EmbedMultipleConfs, kwargs=kwargs)
    except:
        logger.warning(f"rdkit timeout for pose {pose_obj.pose_number}")
        pass
    
    # make obj for protac pose
    protac_pose_obj = classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)
    # for each conformation, align and write to single sdf file,
    # capturing only the poses that obey the rmsd tolerance
    protac_file = linker_folder / 'protac_embedded_confs.sdf'

    try:
        ps = AllChem.MMFFGetMoleculeProperties(protac_embed)
        linker_confs = []

        with open(protac_file, 'a+') as confs_file:
            for i in range(params['rdkit_number_of_confs']):

                # make linker conf "object"
                linker_conf = {'conf_number':i, 'active':None}
                rmsd = Chem.rdMolAlign.AlignMol(protac_embed, reference_ligs_rotate, atomMap=protac_obj.parameters['alignment_pairs'], prbCid=i)
                logger.debug(f"pose: {pose_obj.pose_number}, conf: {i}, rmsd: {rmsd}")

                if rmsd <= params['rmsd_tolerance']:
                    # get internal energy for this conformation
                    ff = AllChem.MMFFGetMoleculeForceField(protac_embed, ps, confId=i)
                    ff.Initialize()
                    linker_conf['energy'] = ff.CalcEnergy()
                    # write this conformation to file
                    if params['write_protac_conf']:
                        molblock = Chem.MolToMolBlock(protac_embed, confId=i, kekulize=False)
                        confs_file.write(f'conf_{i}')
                        confs_file.write(molblock)
                        confs_file.write('$$$$\n')
                    linker_conf['active'] = True
                else:
                    linker_conf['active'] = False
                    linker_conf['energy'] = None
                
                linker_confs.append(linker_conf)

        # if all linker confs were generated successfully, we create their objects:
        for lconf in linker_confs:
            linker_conf = classes.LinkerConf(parent=protac_pose_obj, conf_number=lconf['conf_number'])
            linker_conf.active = lconf['active']
            linker_conf.energy = lconf['energy']
        protac_pose_obj.active = True
        protac_pose_obj.file = protac_file
    except:
        logger.debug(f'No conformation possible for protac {protac_obj.name} on pose {pose_obj.pose_number}')
        protac_pose_obj.active = False
        protac_pose_obj.file = None




