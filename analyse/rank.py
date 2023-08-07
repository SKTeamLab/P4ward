import os
from ..tools.logger import logger
from ..tools import decorators

@decorators.track_run
def protein_poses(
                    ligase_obj,
                    final_ranking_megadock_score,
                    final_ranking_z_score,
                    top_poses,
                    cluster_rep,
                    rank_cluster_reps_only
):
    """
    organize protein pose object in a ranking as specified by the user.
    modifies attibute `active`, which specifies which poses go to
    next stage of development e.g. linker sampling.
    takes all arguments from protein_ranking section in the config file
    """

    rankings = {   
        'final_ranking_megadock_score':{'user_choice':final_ranking_megadock_score, 'ascending':False},
        'final_ranking_z_score':{'user_choice':final_ranking_z_score, 'ascending':True}
    }
    ranking_score_config = [i for i in rankings if rankings[i]['user_choice']][-1]
    ranking = ranking_score_config.replace('final_ranking_','')
    ascending = rankings[ranking_score_config]['ascending']
    # ^ capture the name of the scoring function as it is in the objects' attribute name

    # sort based on the chosen score and save sorted list
    sorted_confs = sorted(ligase_obj.conformations, key=lambda x: getattr(x, ranking), reverse=not ascending)
    ligase_obj.conformations = sorted_confs

    # now sorted list and pose_objs must include only actives
    # pose_objs = [i for i in ligase_obj.conformations if i.active]
    sorted_confs = [i for i in sorted_confs if i.active]
    
    # assign cluster representatives:
    if cluster_rep == 'centroid':
        for i in sorted_confs:
            if i.centroid:
                i.clrep = True
            else: i.clrep = False
    elif cluster_rep == 'best':
        # the best scoring cluster member will appear first
        clusters_ = []
        for i in sorted_confs:
            if i.cluster not in clusters_:
                i.clrep = True
                clusters_.append(i.cluster)
            else:
                i.clrep = False


    if rank_cluster_reps_only:
        # previously the actives came from previous function.
        # now if we must consider only the cluster reps for next stage,
        # then only those should be active
        sorted_confs = [i for i in sorted_confs if i.clrep]
        for pose_obj in ligase_obj.conformations:
            if pose_obj in sorted_confs:
                pose_obj.active = True
            else:
                pose_obj.active = False

    # grab top poses
    sorted_confs = sorted_confs[:top_poses]

    for pose_obj in ligase_obj.conformations:
        if pose_obj not in sorted_confs:
            pose_obj.top = False
        else:
            pose_obj.top = True

    # by the end of the ranking process we have:
    #   - ligase.conformations sorted by chosen score
    #   - then filtered by only the active confs that came from previous function (probably filtering)
    #   - assigned protein cluster representatives based on user choice 'best' or 'centroid'
    #   - if user wants to consider clustering for ranking, deactivated not cluster reps
    #   - from this last version of the rank, grabbed top poses
    #   - assigned the `top` attribute to the top poses. Note that non-top poses should not be deactivated.

    """log user choices"""
    logger.info('generated a final protein-protein ranking:')
    logger.info(f'poses scored by {ranking}')
    if rank_cluster_reps_only:
        logger.info(f'grabbing top {top_poses} poses which are also cluster representatives')
    """"""
    

# @decorators.track_run
def generate_protein_poses(poses, pose_objs, generated_poses_folder, altlocA):
    """
    Subset of protein-protein poses to generate based on conf option `generate_poses`.
    Options are "none", "all", "filtered", "top", "filtered_centroids", "top_centroids".
    Biopython is used to rotate the pose and then save it.
    """
    from Bio.PDB.PDBIO import Select, PDBIO
    pdbio = PDBIO()

    class NotDisordered(Select): # save only altloc A
        def accept_atom(self, atom):
            return not atom.is_disordered() or atom.get_altloc() == "A"

    if poses == "none":
        pass
    elif poses == "all":
        final_poses = pose_objs
    elif poses == 'filtered':
        final_poses = [i for i in pose_objs if i.filtered]
    elif poses == 'top':
        final_poses = [i for i in pose_objs if i.top]
    elif poses == 'filtered_centroids':
        # TODO CHANGE THIS TO REFLECT NEW CLREP
        final_poses = [i for i in pose_objs if i.filtered and i.centroid]
    elif poses == 'top_centroids':
        final_poses = [i for i in pose_objs if i.top and i.centroid]
    else:
        raise Exception("Invalid choice for generate_poses")
    
    generated_poses_folder.mkdir(exist_ok=True)
    for pose in final_poses:
        struct = pose.get_rotated_struct(struct_type='protein', struct_attr='file')

        pdbio.set_structure(struct)
        final_file = generated_poses_folder / f"pose{pose.pose_number}.pdb"
        pose.file = final_file

        if altlocA:
            pdbio.save(str(final_file), select=NotDisordered())
        else:
            pdbio.save(str(final_file))


# @decorators.track_run
def protac_conformations(protac_poses):
    """
    rank the linker conformations for each protac pose and mark the ones that have negative scores.
    """

    for protac_pose in protac_poses:

        def get_linker_score(linker_conf):
            if hasattr(linker_conf, 'rx_score'):
                return(getattr(linker_conf, 'rx_score'))
            else:
                return(None)

        sorted_linker_confs = sorted(protac_pose.linker_confs, key=get_linker_score)
        # sorted_linker_confs = sorted(protac_pose.linker_confs, key=lambda x: getattr(x, 'rx_score'))
        protac_pose.linker_confs = sorted_linker_confs

        print(protac_pose.protein_parent.pose_number, protac_pose.protein_parent.top, [i.rx_score for i in sorted_linker_confs])