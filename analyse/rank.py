from ..tools.logger import logger
from ..tools import decorators

@decorators.track_run
def protein_poses(
                    pose_objs,
                    final_ranking_megadock_score,
                    final_ranking_z_score,
                    top_poses,
                    top_poses_from_centroids_only,
                    use_only_cluster_centroids
):
    """
    organize protein pose object in a ranking as specified by the user.
    modifies attibute `active`, which specifies which poses go to
    next stage of development e.g. linker sampling.
    takes all arguments from protein_ranking section in the config file

    MAYBE, LATER: if more than one scoring function is chosen, the function will
    calculate the consensus between them
    """

    rankings = {   
        'final_ranking_megadock_score':{'user_choice':final_ranking_megadock_score, 'ascending':False},
        'final_ranking_z_score':{'user_choice':final_ranking_z_score, 'ascending':True}
    }
    ranking_score_config = [i for i in rankings if rankings[i]['user_choice']][-1]
    ranking = ranking_score_config.replace('final_ranking_','')
    ascending = rankings[ranking_score_config]['ascending']
    # ^ capture the name of the scoring function as it is in the objects' attribute name

    # sort based on the chosen score
    sorted_confs = sorted(pose_objs, key=lambda x: getattr(x, ranking), reverse=~ascending)
    
    if top_poses_from_centroids_only:
        sorted_confs = [i for i in sorted_confs if i.centroid]

    sorted_confs = sorted_confs[:top_poses]

    if use_only_cluster_centroids:
        sorted_confs = [i for i in sorted_confs if i.centroid]
    
    for pose_obj in pose_objs:
        if pose_obj not in sorted_confs:
            pose_obj.active = False

    """log user choices"""
    logger.info('generated a final protein-protein ranking:')
    logger.info(f'poses scored by {ranking}')
    if top_poses_from_centroids_only:
        logger.info(f'grabbing top {top_poses} poses which are also cluster centroids')
    if use_only_cluster_centroids:
        logger.info(f'grabbing top {top_poses} poses but leaving only the cluster centroids')
    if top_poses_from_centroids_only and use_only_cluster_centroids:
        logger.warning(
            'use_only_cluster_centroids was set as true, but '+
            'top_poses_from_centroids_only supercedes it.'
        )
    """"""
    
