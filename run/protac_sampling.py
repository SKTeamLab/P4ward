from  multiprocessing import Queue, Lock, Process
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
                        linker_scoring_folder,
                        minimize_protac,
                        num_parallel_procs,
                        extend_top_poses_score,
                        extend_top_poses_sampled=False,
                        # pose_objs=None
):

    # information handling:
    # make dict where protac information will be stored and sent to child functions in multiprocs.
    global_parameters = {
        'rdkit_number_of_confs' : rdkit_number_of_confs,
        'protac_poses_folder'   : protac_poses_folder,
        'rmsd_tolerance'        : rmsd_tolerance,
        'time_tolerance'        : time_tolerance,
        'linker_scoring_folder' : linker_scoring_folder,
        'minimize_protac'       : minimize_protac,
        'receptor_obj_file'     : receptor_obj.file
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
    successful_poses = []
    failed_poses = []
    top = len(candidate_poses)

    inQ  = Queue()
    outQ = Queue()
    lock = Lock()

    # send top poses as candidates to inQ
    for pose_obj in candidate_poses:
        params = {
            "pose_number"     : pose_obj.pose_number,
            "pose_obj_rotate" : pose_obj.rotate,
            "file"            : pose_obj.file
        }
        params = {**global_parameters, **params}
        inQ.put(params)
    
    # start processes
    procs = []
    for i in range(num_parallel_procs):
        p = Process(name=i, target=protac_run.sample_protac_pose, args=(inQ, outQ, lock))
        p.daemon = True
        p.start()
        procs.append(p)
    
    # start watching the outQ and counting successful poses
    while len(successful_poses) < top:

        params = outQ.get()

        # rebuild objects from params
        pose_obj = [i for i in ligase_obj.conformations if i.pose_number == params['pose_number']][0]
        protac_pose_obj = classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)
        protac_pose_obj.active = params['protac_pose']['active']
        protac_pose_obj.file   = params['protac_pose']['file']
        try:
            protac_pose_obj.scored_file = params['protac_pose']['scored_file']
        except:
            pass
        
        for linker_conf_dict in params['linker_confs'].values():
            linker_conf = classes.LinkerConf(parent=protac_pose_obj, conf_number=linker_conf_dict['conf_number'])
            linker_conf.active = linker_conf_dict['active']
            linker_conf.rx_score = linker_conf_dict['rx_score']

        success = True
        if extend_top_poses_sampled:
            # if the pose is inactive or has no active linker, success = False, else success = True
            if not params['protac_pose']['active'] or len(params['linker_confs']) == 0:
                success = False
        if extend_top_poses_score:
            # if all the scores of the linker confs are positive, success = False, else success = True
            pos_scores = [True if i['rx_score'] > 0 else False for i in params['linker_confs'].values()]
            if all(pos_scores) or len(pos_scores) == 0:
                success = False

        if success: successful_poses.append(pose_obj)
        else: failed_poses.append(pose_obj)
        
        # get next candidate
        for i in pose_objs:
            if i not in candidate_poses and i not in successful_poses and i not in failed_poses:
                next_candidate = i
                break
        # send it to be processed
        inQ.put(next_candidate)
    
    
    for p in procs:
        p.join()
