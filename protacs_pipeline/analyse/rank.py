from ..tools.logger import logger
from ..tools import decorators

@decorators.track_run
def protein_poses(
                    ligase_obj,
                    top_poses,
):
    """
    organize protein pose object in a ranking as specified by the user.
    modifies attibute `active`, which specifies which poses go to
    next stage of development e.g. linker sampling.
    takes all arguments from protein_ranking section in the config file
    """

    # sort based on the chosen score and save sorted list
    sorted_confs = sorted(ligase_obj.conformations, key=lambda x: getattr(x, 'megadock_score'), reverse=True)
    ligase_obj.conformations = sorted_confs

    # now sorted list and pose_objs must include only actives
    # pose_objs = [i for i in ligase_obj.conformations if i.active]
    sorted_confs = [i for i in sorted_confs if i.active]
    
    # grab top poses
    sorted_confs = sorted_confs[:top_poses]

    for pose_obj in ligase_obj.conformations:
        if pose_obj not in sorted_confs:
            pose_obj.top = False
        else:
            pose_obj.top = True


    """log user choices"""
    logger.info('generated a final protein-protein ranking by megadock score')
    """"""


@decorators.user_choice
@decorators.track_run
def rescore(protac_objs):
    """
    For each protac, rescore its protein poses to combine ppi score
    with protac-protein interaction score (lower=better)
    """

    import numpy as np
    from sklearn.preprocessing import normalize

    for protac_obj in protac_objs:

        protac_pose_objs = protac_obj.active_poses()
        ppis = []
        interactions = []

        for protac_pose_obj in protac_pose_objs:

            interaction = np.mean([i.rx_score for i in protac_pose_obj.active_confs()])
            ppi = protac_pose_obj.protein_parent.megadock_score
            interactions.append(interaction)
            ppis.append(ppi)

        ppis = (np.array(ppis) * -1).reshape(-1,1) 
        # ^ make values negative so that the lower the better to match with interactions
        interactions = np.array(interactions).reshape(-1,1)
        # ^ they are reshaped since they are features so they must be columns

        data = np.concatenate((ppis, interactions), axis=1)
        if len(data) <= 1:
            logger.info('Cannot rescore complexes because only complex was sent to scoring function.')
            return(None)
        data_norm = normalize(data, norm='l2', axis=1)
        scores = np.mean(data_norm, axis=1)

        # record the new scores on the protac pose objs
        for i in range(len(protac_pose_objs)):
            protac_pose_obj = protac_pose_objs[i]
            score = scores[i]
            protac_pose_obj.rescore = score
        
        # sort the protein_poses list for the protac_obj based on the rescores
        def get_rescore(pose_obj):
            return(protac_obj.get_pose(pose_obj).rescore)
        sorted_confs = sorted(protac_obj.protein_poses, key=lambda x: get_rescore(x), reverse=False)
        protac_obj.protein_poses = sorted_confs


@decorators.track_run
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
    else:
        raise Exception("Invalid choice for generate_poses")
    
    generated_poses_folder.mkdir(exist_ok=True)
    for pose in final_poses:
        struct = pose.get_rotated_struct(struct_type='protein')

        pdbio.set_structure(struct)
        final_file = generated_poses_folder / f"pose{pose.pose_number}.pdb"
        pose.file = final_file

        if altlocA:
            pdbio.save(str(final_file), select=NotDisordered())
        else:
            pdbio.save(str(final_file))

