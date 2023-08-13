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

        logger.debug(f"(proc. {p+1} got pose {params['pose_number']})")
        params = protac_run.conf_sampling(params, pose_obj, protac_obj, logger)

        if pose_obj.protac_pose.active:
            print('>>> ', len(pose_obj.protac_pose.linker_confs), pose_obj.protac_pose.linker_confs)

        if pose_obj.protac_pose.active:
            if any(i.active for i in pose_obj.protac_pose.linker_confs):
                params = protac_scoring.rxdock_rescore(params, pose_obj, receptor_obj)
                params = protac_scoring.capture_rxdock_scores(params, pose_obj)
            
            print([i.rx_score for i in pose_obj.protac_pose.active_confs()])
        
        logger.info(f"(proc. {p+1}) Sampled protac for protein pose {params['pose_number']}")
        
        inQ.task_done()
        with lock:
            outQ.put(params)



@decorators.user_choice
# @decorators.track_run
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
                        top_poses,
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
        logger.info(f"Sending pose {params['pose_number']} to protac sampling.")
        inQ.put((params, pose_obj))
    
    # start threads
    threads = []
    for i in range(num_parallel_procs):
        t = Thread(name=i, target=sample_protac_pose, args=(inQ, outQ, lock, i, protac_obj, receptor_obj, logger))
        t.daemon = True
        t.start()
        threads.append(t)
    
    # start watching the outQ and counting successful poses
    while len(successful_poses) <= top_poses:
        print('successful poses:', len(successful_poses))

        params = outQ.get()

        # rebuild objects from params
        # pose_obj = [i for i in ligase_obj.conformations if i.pose_number == params['pose_number']][0]
        # protac_pose_obj = classes.ProtacPose(parent=protac_obj, protein_parent=pose_obj)
        # protac_pose_obj.active = params['protac_pose']['active']
        # protac_pose_obj.file   = params['protac_pose']['file']
        # print('protac_pose active in dict: ', params['protac_pose']['active'])
        # print('protac_pose active in obj: ', protac_pose_obj.active)
        # try:
        #     protac_pose_obj.scored_file = params['protac_pose']['scored_file']
        # except:
        #     pass
        
        # for linker_conf_dict in params['linker_confs'].values():
        #     linker_conf = classes.LinkerConf(parent=protac_pose_obj, conf_number=linker_conf_dict['conf_number'])
        #     linker_conf.active = linker_conf_dict['active']
        #     try:
        #         linker_conf.rx_score = linker_conf_dict['rx_score']
        #     except:
        #         linker_conf.rx_score = None

        # success = True
        success = False
        # if extend_top_poses_sampled:
        #     # if the pose is inactive or has no active linker, success = False, else success = True
        #     if not params['protac_pose']['active'] or len(params['linker_confs']) == 0:
        #         success = False
        # if extend_top_poses_score:
        #     # if all the scores of the linker confs are positive, success = False, else success = True
        #     pos_scores = []
        #     for i in params['linker_confs'].values():
        #         if 'rx_score' in i.keys():
        #             if i['rx_score'] > 0 or i['rx_score'] == None:
        #                 pos_scores.append(True)
        #             else:
        #                 pos_scores.append(False)
        #         else:
        #             pos_scores.append(True)

        #     print(pos_scores)
        #     if all(pos_scores) or len(pos_scores) == 0:
        #         success = False

        if success:
            successful_poses.append(pose_obj)
            pose_obj.top = True
            print('poses: ', [j.conf_number for j in pose_obj.protac_pose.linker_confs])
            print('scores: ', [j.rx_score for j in pose_obj.protac_pose.linker_confs])
        else:
            failed_poses.append(pose_obj)
            pose_obj.top = False

        logger.info(f"pose {params['pose_number']} - {('success' if success else 'failed')}")
        
        # check if all poses have been tried
        if len(candidate_poses) + len(successful_poses) + len(failed_poses) == len(pose_objs):
            logger.info("All possible poses have been sampled.")
            break
            
        # get next candidate
        for i in pose_objs:
            if i not in candidate_poses and i not in successful_poses and i not in failed_poses:
                next_candidate = i
                params = {
                    "pose_number"     : next_candidate.pose_number,
                    "pose_obj_rotate" : next_candidate.rotate,
                    "file"            : next_candidate.file
                }
                params = {**global_parameters, **params}
                break
        
        # send it to be processed
        candidate_poses.append(next_candidate)
        logger.info(f"Sending pose {params['pose_number']} to protac sampling.")
        inQ.put((params, next_candidate))

