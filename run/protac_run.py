import copy
from rdkit.Geometry.rdGeometry import Point3D
from rdkit import Chem
from rdkit.Chem import AllChem
from func_timeout import func_timeout
from ..tools.logger import logger
from ..run.megadock import rotate_atoms


def conf_sampling(params):

    # make folder for each pose obj linker file to be saved
    linker_folder = params['protac_poses_folder'] / f'protein_pose_{params["pose_number"]}'
    linker_folder.mkdir(exist_ok=True)

    # open ligase ligand - changes every loop
    liglig_rotate = copy.deepcopy(params['liglig'])
    conf = liglig_rotate.GetConformer()

    # get final rotation information for the pose
    pose_rotation = params['pose_obj_rotate']

    # rotate ligase ligand to new position
    for i in range(liglig_rotate.GetNumAtoms()):
        x, y, z = conf.GetAtomPosition(i)
        newX, newY, newZ = rotate_atoms((x, y, z), ref_rotation=params['ref_rotation'], pose_rotation=pose_rotation)
        conf.SetAtomPosition(i,Point3D(newX, newY, newZ))

    # combine both ligs
    reference_ligs_rotate = Chem.CombineMols(liglig_rotate, params['reclig'])

    # copy protac to be embedded
    protac_embed = copy.deepcopy(params['protac'])

    # build the dict that will contain the mapping btwn indices and coords:
    coordmap = {}
    mol = reference_ligs_rotate.GetConformer()
    for i in range(len(params['matches']['receptor_lig_indices'])):
        atom_ix = params['matches']['receptor_lig_indices'][i]
        atom_coord = params['matches']['receptor_lig_coords'][i]
        x, y, z = mol.GetAtomPosition(atom_coord)
        coordmap[atom_ix] = Point3D(x, y, z)
    for i in range(len(params['matches']['pose_lig_indices'])):
        atom_ix = params['matches']['pose_lig_indices'][i]
        atom_coord = params['matches']['pose_lig_coords'][i]
        x, y, z = mol.GetAtomPosition(atom_coord)
        coordmap[atom_ix] = Point3D(x, y, z)
    
    # sample protac conformations!
    kwargs = {
        'mol':protac_embed, 'coordMap':coordmap,
        'numConfs':params['rdkit_number_of_confs'],
        'enforceChirality':False
    }
    try:
        func_timeout(params['time_tolerance'], AllChem.EmbedMultipleConfs, kwargs=kwargs)
    except:
        logger.warning(f"rdkit timeout for pose {params['pose_number']}")
        pass

    # for each conformation, align and write to single sdf file,
    # capturing only the poses that obey the rmsd tolerance
    protac_file = linker_folder / 'protac_embedded_confs.sdf'
    params['linker_confs'] = {}
    try:
        with open(protac_file, 'a+') as confs_file:
            for i in range(params['rdkit_number_of_confs']):
                # make linker conf "object"
                linker_conf = {'conf_number':i, 'active':None}
                rmsd = Chem.rdMolAlign.AlignMol(protac_embed, reference_ligs_rotate, atomMap=params['alignment_pairs'], prbCid=i)
                logger.debug(f"pose: {params['pose_number']}, conf: {i}, rmsd: {rmsd}")
                if rmsd <= params['rmsd_tolerance']:
                    molblock = Chem.MolToMolBlock(protac_embed, confId=i, kekulize=False)
                    confs_file.write(f'conf_{i}')
                    confs_file.write(molblock)
                    confs_file.write('$$$$\n')
                    linker_conf['active'] = True
                else:
                    linker_conf['active'] = False
                params['linker_confs'][i] = linker_conf
        params['protac_pose'] = {'active':True, 'file':protac_file}
    except:
        logger.info(f'No conformation possible for pose {params["pose_number"]}')
        params['protac_pose'] = {'active':False, 'file':None}

    # new variables in params dict:
    #   - ['protac_pose']
    #   - ['protac_pose']['file']
    #   - ['protac_pose']['active']
    #   - ['linker_confs']

    return(params)



def sample_protac_pose(inQ, outQ, lock, p):

    from . import protac_scoring

    while True:

        with lock:
            params = inQ.get()
            logger.debug(f"(proc. {p+1} got pose {params['pose_number']})")
        
        params = conf_sampling(params)

        if (
            params['protac_pose']['active'] and  # only rescore if protac pose is active
            len(params['linker_confs']) > 0 and  # if it generated at least one linker conf
            any([i['active'] for i in params['linker_confs'].values()]) # if at least one is active
        ):
            params = protac_scoring.rxdock_rescore(params)
            params = protac_scoring.capture_rxdock_scores(params)
        
        logger.info(f"(proc. {p+1}) Sampled protac for protein pose {params['pose_number']}")
        
        with lock:
            outQ.put(params)
