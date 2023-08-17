from  threading import Lock, Thread
from queue import Queue
from ..tools.logger import logger
from ..tools import classes
from ..tools import decorators


def sample_protac_pose(inQ, outQ, lock, p, protac_obj, receptor_obj, logger):

    from . import protac_scoring
    from . import protac_run

    while True:

        with lock:
            params, pose_obj = inQ.get()

        logger.debug(f"(proc. {p+1} got pose {pose_obj.pose_number})")
        protac_run.conf_sampling(params, pose_obj, protac_obj, logger)

        if pose_obj.protac_pose.active:
            if any(i.active for i in pose_obj.protac_pose.linker_confs):
                protac_scoring.rxdock_rescore(params, pose_obj, receptor_obj)
                protac_scoring.capture_rxdock_scores(pose_obj)
            
        logger.info(f"(proc. {p+1}) Sampled protac for protein pose {pose_obj.pose_number}")
        
        inQ.task_done()
        with lock:
            outQ.put(pose_obj)



@decorators.user_choice
@decorators.track_run
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
        'minimize_protac'       : minimize_protac
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

    # ~~~~~~~~~~~~~~
    # start sampling
    # ~~~~~~~~~~~~~~

    from . import protac_run

    pose_objs = [i for i in ligase_obj.conformations if i.active]
    candidate_poses = [i for i in pose_objs if i.top]
    top_poses = len(candidate_poses)
    successful_poses = []
    failed_poses = []

    inQ  = Queue()
    outQ = Queue()
    lock = Lock()

    # send top poses as candidates to inQ
    for pose_obj in candidate_poses:
        logger.info(f"Sending pose {pose_obj.pose_number} to protac sampling.")
        inQ.put((global_parameters, pose_obj))
    
    # start threads
    threads = []
    for i in range(num_parallel_procs):
        t = Thread(name=i, target=sample_protac_pose, args=(inQ, outQ, lock, i, protac_obj, receptor_obj, logger))
        t.daemon = True
        t.start()
        threads.append(t)
    
    # start watching the outQ and counting successful poses
    while len(successful_poses) <= top_poses:

        pose_obj = outQ.get()  
        # sort linker conformations based on score
        pose_obj.protac_pose.linker_confs = sorted(pose_obj.protac_pose.active_confs(), key=lambda x: getattr(x, 'rx_score'))

        success = True
        if extend_top_poses_sampled:
            # if the pose is inactive or has no active linker, success = False, else success = True
            if pose_obj.protac_pose.active == False or len(pose_obj.protac_pose.active_confs()) == 0:
                success = False
            else:
                if extend_top_poses_score:
                    # if all the scores of the linker confs are positive, success = False, else success = True
                    pos_scores = []
                    for i in pose_obj.protac_pose.active_confs():
                        if i.rx_score > 0:
                            pos_scores.append(True)
                        else:
                            pos_scores.append(False)
                    if all(pos_scores) or len(pos_scores) == 0:
                        success = False
                else:
                    success=True


        if success:
            successful_poses.append(pose_obj)
            pose_obj.top = True
        else:
            failed_poses.append(pose_obj)
            pose_obj.top = False

        logger.info(f"pose {pose_obj.pose_number} - {('success' if success else 'failed')}")
        
        # check if all poses have been tried
        if len(candidate_poses) + len(successful_poses) + len(failed_poses) == len(pose_objs):
            logger.info("All possible poses have been sampled.")
            break
            
        # get next candidate
        for i in pose_objs:
            if i not in candidate_poses and i not in successful_poses and i not in failed_poses:
                next_candidate = i
                break
        
        # send it to be processed
        candidate_poses.append(next_candidate)
        logger.info(f"Sending pose {next_candidate.pose_number} to protac sampling.")
        inQ.put((global_parameters, next_candidate))

