import os
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
    """

    active_pose_objs = [i for i in pose_objs if i.active]

    rankings = {   
        'final_ranking_megadock_score':{'user_choice':final_ranking_megadock_score, 'ascending':False},
        'final_ranking_z_score':{'user_choice':final_ranking_z_score, 'ascending':True}
    }
    ranking_score_config = [i for i in rankings if rankings[i]['user_choice']][-1]
    ranking = ranking_score_config.replace('final_ranking_','')
    ascending = rankings[ranking_score_config]['ascending']
    # ^ capture the name of the scoring function as it is in the objects' attribute name

    # sort based on the chosen score
    sorted_confs = sorted(active_pose_objs, key=lambda x: getattr(x, ranking), reverse=~ascending)
    
    if top_poses_from_centroids_only:
        sorted_confs = [i for i in sorted_confs if i.centroid]

    sorted_confs = sorted_confs[:top_poses]

    if use_only_cluster_centroids:
        sorted_confs = [i for i in sorted_confs if i.centroid]
    
    for pose_obj in pose_objs:
        if pose_obj not in sorted_confs:
            pose_obj.active = False
            pose_obj.top = False
        else:
            pose_obj.top = True

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
    


def generate_protein_poses(poses, pose_objs, generated_poses_folder):
    """
    Subset of protein-protein poses to generate based on conf option `generate_poses`.
    Options are "none", "all", "filtered", "top", "filtered_centroids", "top_centroids".
    Biopython is used to rotate the pose and then save it.
    """
    from Bio.PDB.PDBIO import PDBIO
    pdbio = PDBIO()

    if poses == "none":
        pass
    elif poses == "all":
        final_poses = pose_objs
    elif poses == 'filtered':
        final_poses = [i for i in pose_objs if i.filtered]
    elif poses == 'top':
        final_poses = [i for i in pose_objs if i.top]
    elif poses == 'filtered_centroids':
        final_poses = [i for i in pose_objs if i.filtered and i.centroid]
    elif poses == 'top_centroids':
        final_poses = [i for i in pose_objs if i.top and i.centroid]
    else:
        raise Exception("Invalid choice for generate_poses")
    
    
    for pose in final_poses:
        struct = pose.get_rotated_struct(struct_type='protein', struct_attr='mg_file')
        
        pdbio.set_structure(struct)
        final_file = os.path.join(generated_poses_folder, f"pose{pose.pose_number}.pdb")
        pdbio.save(final_file)
