import subprocess
import os
from ..tools import decorators
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

    pymol_combine(receptor_file, receptor_ligand_file, out_filename=prep_receptor_file)
    pymol_combine(ligase_file,   ligase_ligand_file,   out_filename=prep_ligase_file)
 
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
    if num_threads == 'all': pass 
    else: os.putenv('OMP_THREAD_LIMIT', num_threads)

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
    generate the ProteinPose objs from the docking run
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
    output.close()


def rotate_atoms(atom_coords, ref_rotation, pose_rotation):
    """
    move atoms as performed by megadock decoygen, using information
    captured from the megadock output file
    """

    ref_psi, ref_theta, ref_phi = ref_rotation['initial_rot']
    l1, l2, l3 = ref_rotation['lig_transl']
    r1, r2, r3 = ref_rotation['rec_transl']
    N = ref_rotation['grid_size']
    spacing = ref_rotation['spacing']

    t1, t2, t3 = pose_rotation['transl']
    a1, a2, a3 = pose_rotation['angles']

    def rotate(psi, theta, phi, oldX, oldY, oldZ):
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

        newX = r11 * oldX + r12 * oldY + r13 * oldZ
        newY = r21 * oldX + r22 * oldY + r23 * oldZ
        newZ = r31 * oldX + r32 * oldY + r33 * oldZ

        return(newX, newY, newZ)
    
    # first subtract ligand initial translation
    coord1, coord2, coord3 = atom_coords
    oldx = coord1 - l1
    oldy = coord2 - l2
    oldz = coord3 - l3

    # rotate based on reference rotation
    tx1, ty1, tz1 = rotate(ref_psi, ref_theta, ref_phi, oldx, oldy, oldz)
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
def cluster(pose_objects, clustering_cutoff):
    """
    cluster megadock docked poses using a user-specified RMSD cutoff
    RMSD is calculated using the alpha carbons.
    """
    # NOTE list is the active ligase.conformations

    from sklearn.cluster import AgglomerativeClustering
    from sklearn.neighbors import NearestCentroid
    from sklearn.neighbors import NearestNeighbors
    from .structure_tools import get_rmsd
    import numpy as np

    logger.info(f'Clustering protein poses using cutoff of {clustering_cutoff}')
    
    # get first pose as reference:
    reference_obj = pose_objects[0]
    reference_obj_struct = reference_obj.get_rotated_struct(struct_type='protein')

    # load receptor object
    for pose_obj in pose_objects:
        pose_obj_struct = pose_obj.get_rotated_struct(struct_type='protein')

        rmsd = get_rmsd(reference_obj_struct, pose_obj_struct, ca=True)
        pose_obj.rmsd = rmsd
        pose_obj.rmsd_reference = reference_obj
    
    # get the array of rmsds to compute
    rmsd_list = [i.rmsd for i in pose_objects]
    rmsd_points = np.asarray(rmsd_list).reshape(-1, 1)

    # run clustering
    clustering = AgglomerativeClustering(
        n_clusters=None,
        linkage='average',
        distance_threshold=clustering_cutoff
    ).fit(rmsd_points)

    # add the information to the objects as attributes
    for i in range(len(pose_objects)):
        pose_objects[i].cluster = clustering.labels_[i]

    # get the theoretical centroid of each cluster
    nc = NearestCentroid()
    nc.fit(rmsd_points, clustering.labels_)
    # get which rmsd value is the centroid
    # the centroid will be the one closest the theoretical centroid
    knn = NearestNeighbors(n_neighbors=1)
    knn.fit(rmsd_points)
    _, indices = knn.kneighbors(nc.centroids_)
    centroids = np.hstack(rmsd_points)[np.hstack(indices)]

    # add to the attributes whether the pose is a centroid or not
    for pose_obj in pose_objects:
        if pose_obj.rmsd.round(decimals=3) in centroids.round(decimals=3):
            pose_obj.centroid = True
        else:
            pose_obj.centroid = False