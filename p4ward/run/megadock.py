import subprocess
import os
from ..tools import decorators, classes
from ..tools.logger import logger
from ..definitions import CWD


@decorators.user_choice # NOTE choice is True if run megadock is true
@decorators.track_run
def prep_structures(receptor_obj, ligase_obj):
    """
    Combine and prepare the protein files with their ligands for docking
    using pymol combine
    """

    from .structure_tools import pymol_combine

    receptor_file = receptor_obj.active_file
    ligase_file = ligase_obj.active_file
    receptor_ligand_file = receptor_obj.lig_file
    ligase_ligand_file = ligase_obj.lig_file

    prep_receptor_file = receptor_file.parent/('mg-'+receptor_file.name)
    prep_ligase_file = ligase_file.parent/('mg-'+ligase_file.name)

    pymol_combine(receptor_file, receptor_ligand_file, out_filename=prep_receptor_file, assign_chains=False)
    pymol_combine(ligase_file,   ligase_ligand_file,   out_filename=prep_ligase_file,   assign_chains=False)
 
    logger.info(f'Saved protein files for megadock: {prep_receptor_file}, {prep_ligase_file}')

    receptor_obj.mg_file = prep_receptor_file
    ligase_obj.mg_file = prep_ligase_file


@decorators.user_choice
@decorators.track_run
def run_docking(program_path, receptor_file, ligase_file, num_threads, run_docking_output_file,
                num_predictions, num_predictions_per_rotation, num_rotational_angles, log_file):
    """
    Run the megadock main step
    """

    # if num_threads not specified, megadock uses all
    if num_threads != 'all':
        os.putenv('OMP_THREAD_LIMIT', num_threads)

    command = [
        program_path,
        '-R', receptor_file,
        '-L', ligase_file,
        '-N', num_predictions,
        '-t', num_predictions_per_rotation,
        '-r', num_rotational_angles,
        '-o', run_docking_output_file
    ]

    logger.info('Running megadock...')
    with open(log_file, 'w+') as log_file_: 
        subprocess.run(command, stdout=log_file_)
    logger.info('Done.')


@decorators.user_choice # NOTE choice is True if run megadock is true
@decorators.track_run
def capture_scores(run_docking_output_file, ligase_obj):
    """
    Grab all original megadock scores from run_docking_output_file and
    generate the ProteinPose objs from the docking run.

    :param run_docking_output_file: file generated my megadock output
    :param ligase_obj: object representing the ligase
    """
    from ..tools.classes import ProteinPose
    logger.info(f'Capturing megadock scores from {run_docking_output_file}')

    # start parsing the output
    ## grab header rotation information and add to lig_obj
    output = open(run_docking_output_file, 'r')
    header = [next(output) for _ in range(4)]
    rotate = {
        'grid_size':   int(header[0].split('\t')[0]),
        'spacing':     float(header[0].split('\t')[1]),
        'initial_rot': [float(i) for i in header[1].split('\t') if i!=''],
        'rec_transl':  [float(i) for i in header[2].split('\t')[-3:]],
        'lig_transl':  [float(i) for i in header[3].split('\t')[-3:]]
    }
    ligase_obj.rotate = rotate
    
    ## grab the information for each pose
    count = 0
    for line in output:
        if len(line.split('\t')) == 7:
            count += 1

            rotate = {
                'angles':[float(i) for i in line.split('\t')[:3]],
                'transl':[float(i) for i in line.split('\t')[3:6]]
            }
            score = float(line.split('\t')[-1])

            pose_obj = ProteinPose(parent=ligase_obj, pose_number=int(count))
            pose_obj.megadock_score = score
            pose_obj.rotate = rotate
            pose_obj.active = True
            pose_obj.filter_info = {}
    output.close()


def rotate_atoms(atom_coords, ref_rotation, pose_rotation):
    """
    Move atoms using information captured from the megadock output file.
    """

    ref_psi, ref_theta, ref_phi = ref_rotation['initial_rot']
    l1, l2, l3 = ref_rotation['lig_transl']
    r1, r2, r3 = ref_rotation['rec_transl']
    N = ref_rotation['grid_size']
    spacing = ref_rotation['spacing']

    t1, t2, t3 = pose_rotation['transl']
    a1, a2, a3 = pose_rotation['angles']

    def rotate(psi, theta, phi, prevX, prevY, prevZ):
        import numpy as np

        r11 = np.cos(psi)*np.cos(phi)  -  np.sin(psi)*np.cos(theta)*np.sin(phi)
        r21 = np.sin(psi)*np.cos(phi)  +  np.cos(psi)*np.cos(theta)*np.sin(phi)
        r31 = np.sin(theta)*np.sin(phi)

        r12 = -np.cos(psi)*np.sin(phi)  -  np.sin(psi)*np.cos(theta)*np.cos(phi)
        r22 = -np.sin(psi)*np.sin(phi)  +  np.cos(psi)*np.cos(theta)*np.cos(phi)
        r32 = np.sin(theta)*np.cos(phi)

        r13 = np.sin(psi)*np.sin(theta)
        r23 = -np.cos(psi)*np.sin(theta)
        r33 = np.cos(theta)

        newX = r11 * prevX + r12 * prevY + r13 * prevZ
        newY = r21 * prevX + r22 * prevY + r23 * prevZ
        newZ = r31 * prevX + r32 * prevY + r33 * prevZ

        return(newX, newY, newZ)
    
    # first subtract ligand initial translation
    coord1, coord2, coord3 = atom_coords
    prevx = coord1 - l1
    prevy = coord2 - l2
    prevz = coord3 - l3

    # rotate based on reference rotation
    tx1, ty1, tz1 = rotate(ref_psi, ref_theta, ref_phi, prevx, prevy, prevz)
    # rotate based on pose angles
    tx2, ty2, tz2 = rotate(a1, a2, a3, tx1, ty1, tz1)

    if t1 >= N/2: t1 -= N
    if t2 >= N/2: t2 -= N
    if t3 >= N/2: t3 -= N

    # calculate to find final box
    final_x = tx2-t1*spacing+r1
    final_y = ty2-t2*spacing+r2
    final_z = tz2-t3*spacing+r3

    return(final_x, final_y, final_z)


@decorators.user_choice
@decorators.track_run
def cluster(ligase_obj, protac_objs, clustering_type, clustering_cutoff_redund, clustering_cutoff_trend, cluster_redund_repr, rescore_poses):
    """
    cluster protein_poses' triad points for redundancy or trend using hieragglo with ward linkage
    """

    import numpy as np
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.neighbors import NearestCentroid
    from sklearn.neighbors import NearestNeighbors


    def make_clusterer(pose_objs, rescore_poses, cutoff):

        from .structure_tools import get_coords_array

        logger.info(f'Clustering protein poses for {clustering_type} using cutoff of {cutoff}')

        ## start clustering
        coords = get_coords_array(pose_objs, ligase_obj)
        
        clusterer = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=cutoff,
            metric='euclidean',
            linkage='ward'
        ).fit(coords)

        ## get the centroids
        if clusterer.n_clusters_ == 1:
            # cannot use NearestCentroid if there is only 1 cluster
            centroids = coords.mean(axis=0).reshape(1,-1)
        else:
            nc = NearestCentroid()
            nc.fit(coords, clusterer.labels_)
            centroids = nc.centroids_

        knn = NearestNeighbors(n_neighbors=1)
        knn.fit(coords)
        _, indices = knn.kneighbors(centroids)

        repr_centr = list(np.asarray(pose_objs)[indices].ravel())
        repr_best = []

        ## make cluster information and get best scores for each cluster
        cluster = classes.Cluster(clusterer=clusterer, type=clustering_type)
        
        for cln in range(clusterer.n_clusters_):
            cl_components = list(np.asarray(pose_objs)[clusterer.labels_ == cln])
            cluster.clusters[cln] = cl_components
            if rescore_poses:
                protac_components = [protac_obj.get_pose(i) for i in cl_components]
                best = protac_components[np.argmin(i.rescore for i in protac_components)]
            else:
                best = cl_components[np.argmax([i.megadock_score for i in cl_components])]
            repr_best.append(best)

        cluster.repr_centr = repr_centr
        cluster.repr_best = repr_best
        cluster.coords = coords
        ## get best scored repr

        return(cluster)

    def enough_poses(pose_objs):
        check = len(pose_objs) > 1
        if not check:
            logger.info("Cannot cluster protein poses because fewer than 2 poses were sent to the clustering step. Skipping.")
        return(check)


    if clustering_type == 'redundancy':
        
        pose_objs = ligase_obj.active_confs()
        if enough_poses(pose_objs):

            cluster_obj = make_clusterer(pose_objs=pose_objs, rescore_poses=rescore_poses, cutoff=clustering_cutoff_redund)
            ligase_obj.cluster = cluster_obj

            if cluster_redund_repr == 'centroid':
                attribute = 'repr_centr'
            elif cluster_redund_repr == 'best':
                attribute = 'repr_best'

            for pose_obj in pose_objs:
                if pose_obj in getattr(cluster_obj, attribute):
                    pose_obj.active = True
                else:
                    pose_obj.active = False
        

    elif clustering_type == 'trend':

        for protac_obj in protac_objs:

            pose_objs = protac_obj.protein_poses
            if enough_poses(pose_objs):
                cluster_obj = make_clusterer(pose_objs=pose_objs, rescore_poses=rescore_poses, cutoff=clustering_cutoff_trend)
                protac_obj.cluster = cluster_obj
