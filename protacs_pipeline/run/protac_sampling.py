from  threading import Thread
from queue import Queue
from ..tools.logger import logger
from ..tools import decorators


def sample_protac_pose(inQ, outQ, p, receptor_obj, ligase_obj, global_parameters, logger):

    from . import protac_scoring
    from . import protac_run

    while True:

        protac_obj, pose_obj = inQ.get()

        logger.debug(f"(proc. {p+1} got pose {pose_obj.pose_number})")
        protac_run.conf_sampling(global_parameters, pose_obj, protac_obj, logger)

        protac_pose = protac_obj.get_pose(pose_obj)
        
        if protac_pose.active == True and global_parameters['rxdock_score']:
            if any(i.active for i in protac_pose.linker_confs):
                protac_scoring.rxdock_rescore(global_parameters, pose_obj, receptor_obj, ligase_obj, protac_obj)
                protac_scoring.capture_rxdock_scores(pose_obj, protac_obj)
            
        logger.debug(f"(proc. {p+1}) Sampled protac {protac_obj.name} for protein pose {pose_obj.pose_number}")
        
        inQ.task_done()
        outQ.put((pose_obj, protac_obj))



@decorators.user_choice
@decorators.track_run
def protac_sampling(
                        receptor_obj,
                        ligase_obj,
                        protac_objs,
                        extend_flexible_small_linker,
                        neighbour_number,
                        min_linker_length,
                        rdkit_number_of_confs,
                        rdkit_ligands_cleanup,
                        write_protac_conf,
                        protac_poses_folder,
                        rmsd_tolerance,
                        time_tolerance,
                        unbound_protac_num_confs,
                        rxdock_score,
                        linker_scoring_folder,
                        minimize_protac,
                        num_parallel_procs,
                        extend_top_poses_score,
                        extend_top_poses_sampled=False,
                        extend_top_poses_energy=False
                        # pose_objs=None
):

    # information handling:
    # make dict where protac information will be stored and sent to child functions in multiprocs.
    global_parameters = {
        'rdkit_number_of_confs' : rdkit_number_of_confs,
        'protac_poses_folder'   : protac_poses_folder,
        'rmsd_tolerance'        : rmsd_tolerance,
        'time_tolerance'        : time_tolerance,
        'rxdock_score'          : rxdock_score,
        'linker_scoring_folder' : linker_scoring_folder,
        'minimize_protac'       : minimize_protac,
        'write_protac_conf'     : write_protac_conf
    }

    # make folder where the linkers for all pose objs will be stored
    protac_poses_folder.mkdir(exist_ok=True)

    # ~~~~~~~~~~~~~~
    # prepare protac
    # ~~~~~~~~~~~~~~

    from .protac_prep import protac_prep

    # if protac unbound energy will be needed but it has not been sampled yet:
    for protac_obj in protac_objs:
        if extend_top_poses_energy and not hasattr(protac_obj, 'unbound_energy'):
            protac_obj.sample_unbound_confs(num_unbound_confs=unbound_protac_num_confs)

    for protac_obj in protac_objs:
        protac_obj.parameters = protac_prep(
            receptor_obj, ligase_obj, protac_obj,
            extend_flexible_small_linker,
            min_linker_length,
            neighbour_number,
            rdkit_ligands_cleanup
        )

    # ~~~~~~~~~~~~~~
    # start sampling
    # ~~~~~~~~~~~~~~

    pose_objs = [i for i in ligase_obj.conformations if i.active]

    inQ  = Queue()
    outQ = Queue()

    # start threads
    threads = []
    for i in range(num_parallel_procs):
        t = Thread(name=i, target=sample_protac_pose, args=(inQ, outQ, i, receptor_obj, ligase_obj, global_parameters, logger))
        t.daemon = True
        t.start()
        threads.append(t)


    for protac_obj in protac_objs:

        logger.info(f"Now sampling protac {protac_obj.name}")

        candidate_poses = [i for i in pose_objs if i.top]
        top_poses = len(candidate_poses)
        successful_poses = []
        failed_poses = []
        
        # send top poses as candidates to inQ
        for pose_obj in candidate_poses:
            logger.debug(f"Sending pose {pose_obj.pose_number} to sample protac {protac_obj.name}.")
            inQ.put((protac_obj, pose_obj))
    
        # start watching the outQ and counting successful poses
        while True:

            logger.info(f"Sampled {len(successful_poses)}/{top_poses} successful poses.")
            if len(successful_poses) == top_poses:
                break

            pose_obj, protac_obj = outQ.get()  

            protac_pose = protac_obj.get_pose(pose_obj)
            # sort linker conformations based on score
            if rxdock_score:
                protac_pose.linker_confs = sorted(protac_pose.active_confs(), key=lambda x: getattr(x, 'rx_score'))

            success = True
            if extend_top_poses_sampled:
                # if the pose is inactive or has no active linker, success = False, else success = True
                if protac_pose.active == False or len(protac_pose.active_confs()) == 0:
                    success = False
                else:
                    if extend_top_poses_score:
                        # if all the scores of the linker confs are positive, success = False, else success = True
                        pos_scores = []
                        for i in protac_pose.active_confs():
                            if i.rx_score > 0:
                                pos_scores.append(True)
                            else:
                                pos_scores.append(False)
                        if all(pos_scores) or len(pos_scores) == 0:
                            success = False
                    else:
                        success=True
                    if extend_top_poses_energy:
                        higher_energies = []
                        for i in protac_pose.active_confs():
                            if i.energy > protac_obj.unbound_energy:
                                higher_energies.append(True)
                            else:
                                higher_energies.append(False)
                        if all(higher_energies) or len(higher_energies) == 0:
                            success = False
                    else:
                        success=True

            if success:
                successful_poses.append(pose_obj)
                protac_obj.protein_poses.append(pose_obj)
            else:
                failed_poses.append(pose_obj)

            logger.debug(f"pose {pose_obj.pose_number}, protac {protac_obj.name} - {('success' if success else 'failed')}")
            
            # check if all poses have been tried
            if len(successful_poses) + len(failed_poses) == len(pose_objs):
                logger.info(f"All possible poses for protac {protac_obj.name} have been sampled.")
                break
            
            # IMPORTANT! we can only send a next candidate if the current one was NOT successful
            # and only if there is a next candidate
            if not success:
                # get next candidate
                for i in pose_objs:
                    if i not in candidate_poses and i not in successful_poses and i not in failed_poses:
                        next_candidate = i
                        break
                    else:
                        next_candidate = None
                
                # send it to be processed
                if next_candidate != None:
                    candidate_poses.append(next_candidate)
                    logger.debug(f"Sending pose {next_candidate.pose_number} to sample protac {protac_obj.name}.")
                    inQ.put((protac_obj, next_candidate))

            outQ.task_done()
                
        logger.info(f"Finished sampling for protac {protac_obj.name}")
    logger.info(f"Sampling done for all {len(protac_objs)} protac linker variations.")