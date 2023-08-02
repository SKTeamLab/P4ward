from  multiprocessing import Queue, Lock
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Geometry.rdGeometry import Point3D
from ..tools.logger import logger
from ..tools import classes
from . import protac_prep


# rdkit_protac_sampling() (here)             = coordinate_everything()
#   make_lig_indices()           |
#   match_substructures()        |-> .protac_prep.protac_prep()
#   check_linker_size()          |
#   
#   for i in protac_pose
#       sample_protac_pose() (here)          = run_steps_for_one_obj()
#           conf_sampling()         |
#           conf_scoring()          |-> .protac_run
#           capture_conf_score()    | 



def protac_sampling(
                        receptor_obj,
                        ligase_obj,
                        protac_obj,
                        extend_flexible_small_linker,
                        neighbour_number,
                        min_linker_length,
                        rdkit_number_of_confs,
                        protac_poses_folder,
                        rmsd_tolerance,
                        time_tolerance,
                        extend_top_poses_sampled=False,
                        # pose_objs=None
):

    # information handling:
    # make dict where protac information will be stored and sent to child functions in multiprocs.
    global_parameters = {
        'rdkit_number_of_confs' : rdkit_number_of_confs,
        'protac_poses_folder'   : protac_poses_folder,
        'rmsd_tolerance'        : rmsd_tolerance,
        'time_tolerance'        : time_tolerance
    }

    # make folder where the linkers for all pose objs will be stored
    protac_poses_folder.mkdir(exist_ok=True)

    # ~~~~~~~~~~~~~~
    # prepare protac
    # ~~~~~~~~~~~~~~

    from .protac_prep import protac_prep
    parameters = protac_prep(
        receptor_obj, ligase_obj, protac_obj,
        extend_flexible_small_linker,
        min_linker_length,
        neighbour_number
    )

    global_parameters.update(parameters)

    # global parameters will later be joined with pose specific parameters
    # and sent to sample_protac_pose()

    # ~~~~~~~~~~~~~~
    # start sampling
    # ~~~~~~~~~~~~~~

    from . import protac_run

    pose_objs = [i for i in ligase_obj.conformations if i.active]
    candidate_poses = [i for i in pose_objs if i.top]
    top = len(candidate_poses)

    q = Queue()

    for pose_obj in candidate_poses:
        params = {
            "pose_number" : pose_obj.pose_number,
            "pose_obj_rotate" : pose_obj.rotate
        }
        params = {**global_parameters, **params}
        q.put(params)
    

    while q.qsize() > 0:
        
        # run all functions for a single protein pose
        params = protac_run.sample_protac_pose(q=q)

        # build objects and their attributes from the returned params

        pose_obj = [i for i in ligase_obj.conformations if i.pose_number == params['pose_number']][0]
        
        protac_pose_obj = classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)
        protac_pose_obj.active = params['protac_pose']['active']
        protac_pose_obj.file   = params['protac_pose']['file']
        
        for linker_conf_dict in params['linker_confs']:
            linker_conf = classes.LinkerConf(parent=protac_pose_obj, conf_number=linker_conf_dict['conf_number'])
            linker_conf.active = linker_conf_dict['active']

        print(f"done with protein_pose {pose_obj.pose_number}")
