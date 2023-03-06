import subprocess
import os
from ..tools import decorators
from ..tools.logger import logger
from ..definitions import CWD


@decorators.user_choice
@decorators.track_run
def run_docking(program_path, receptor_file, ligase_file, num_threads, run_docking_output_file,
                num_predictions, num_predictions_per_rotation, log_file):
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
        '-o', run_docking_output_file
    ]

    logger.info('Running megadock...')
    with open(log_file, 'w+') as log_file_: 
        subprocess.run(command, stdout=log_file_)
    logger.info('Done.')


@decorators.user_choice # choice is True if run megadock is true
def capture_scores(run_docking_output_file, ligase_obj):
    """
    Grab all original megadock scores from run_docking_output_file and
    generate the ProteinPose objs from the docking run
    """
    from ..tools.classes import ProteinPose
    logger.info(f'Capturing megadock scores from {run_docking_output_file}')

    # start parsing the output
    output = open(run_docking_output_file, 'r')
    count = 0
    for line in output:
        if len(line.split('\t')) == 7:

            count += 1
            score = float(line.split('\t')[-1])
            pose_obj = ProteinPose(parent=ligase_obj, pose_number=int(count))
            pose_obj.megadock_score = score
            pose_obj.active = True
    output.close()


@decorators.user_choice
@decorators.track_run
def generate_poses(run_docking_output_file, ligase_obj, docked_poses_folder):
    """
    If user chooses to generate all poses from megadock, independent of 
    filtering, this function does it and adds the file attribute.
    """
    from ..tools.script_tools import create_folder

    logger.info(
        f'Generating all poses from output file {run_docking_output_file} '+
        f'and saving them in the folder {docked_poses_folder}.'
    )

    create_folder(docked_poses_folder)

    output = open(run_docking_output_file, 'r')
    count = 0
    for line in output:
        if len(line.split('\t')) == 7:

            count += 1
            pose_obj = ligase_obj.conformations[count - 1]

            decoy_name = os.path.join(docked_poses_folder, f'decoy{count}.pdb')
            decoygen_command = [
                'decoygen', decoy_name,
                ligase_obj.file,
                run_docking_output_file, str(count)
            ]
            subprocess.run(decoygen_command, stdout=subprocess.DEVNULL)
            pose_obj.file = decoy_name

            # check if file was created successfully
            if not os.path.isfile(pose_obj.file):
                raise Exception("Could not generate pose successfully")

    output.close()


@decorators.user_choice
@decorators.track_run
def filter_poses(receptor_obj, ligase_obj, dist_cutoff,
                 output_file, output_filtered_file, docked_poses_folder,
                 generate_all_poses=False):
    """
    Use Biopython to filter the megadock poses which satisfy a
    distance cutoff for both binding sites. Takes the a Protein object
    and handles the rest by accessing its attributes.
    """

    from ..tools.structure_tools import structure_proximity
    from ..tools.script_tools import create_folder

    logger.info(f'Filtering megadock poses with cuttoff {dist_cutoff}')

    output = open(output_file, 'r')
    output_filtered = open(output_filtered_file, 'a+')

    # make a folder for the docked structures
    create_folder(docked_poses_folder)

    count = 0
    for line in output:

        if len(line.split('\t')) == 7:
            count += 1
            ligase_pose_obj = ligase_obj.conformations[count - 1]

            # we will only generate the poses that were not already generated before
            if ligase_pose_obj.file is None or not os.path.isfile(ligase_pose_obj.file):
                decoy_name = os.path.join(docked_poses_folder, f'decoy{count}.pdb')
                decoygen_command = [
                    'decoygen', decoy_name,
                    ligase_obj.file,
                    output_file, str(count)
                ]
                subprocess.run(decoygen_command, stdout=subprocess.DEVNULL)
                ligase_pose_obj.file = decoy_name

            distance, proximity = structure_proximity(
                receptor_obj.ligand_struct,
                ligase_pose_obj.get_ligand_struct(),
                dist_cutoff=dist_cutoff
            )

            if proximity:
                output_filtered.write(line)
                ligase_pose_obj.active = True
                logger.info(f'Activating filtered pose {count}, with distance of {distance} between ligands.')
            else:
                ligase_pose_obj.active = False
                if not generate_all_poses:
                    ligase_pose_obj.file = None
                    os.remove(decoy_name)

        else: 
            output_filtered.write(line)
    
    logger.info(f'Saved filtered megadock poses to file {output_filtered_file}')

    output.close()
    output_filtered.close()


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
    from ..tools.structure_tools import get_rmsd
    import numpy as np
    
    # get first pose as reference:
    reference_obj = pose_objects[0]
    reference_obj_struct = reference_obj.get_protein_struct()

    # load receptor object
    for pose_obj in pose_objects:
        pose_obj_struct = pose_obj.get_protein_struct()

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


@decorators.user_choice
@decorators.track_run
def zrank_rescore(ligase_obj, receptor_obj, zrank_path, run_docking_output_file):
    """
    use zrank to rescore a previous "pure" megadock run
    """
    import tempfile
    from ..tools.structure_tools import reduce

    # check if receptor and ligase were already reduced
    for protein_obj in (ligase_obj, receptor_obj):
        if not hasattr(protein_obj, 'reduced_file'):
            reduce([protein_obj])

    # zrank considers all conformations + 1
    conf_count = len(ligase_obj.conformations)
    # read megadock outpu
    megadock_output_read = open(run_docking_output_file, 'r').read()
    # the end of the file has an empty line, zrank tries to score that. So strip it out:
    megadock_output_read = megadock_output_read.rstrip()
    # now swap protein file names to reduced file names
    megadock_output_read = megadock_output_read.replace(receptor_obj.file, receptor_obj.reduced_file)
    megadock_output_read = megadock_output_read.replace(ligase_obj.file, ligase_obj.reduced_file)
    
    # run zrank using `megadock_output_read` as temporary file
    with tempfile.NamedTemporaryFile(mode='w', dir=CWD, delete=True) as tmp:
        tmp.write(megadock_output_read)
        command = [
            zrank_path,
            tmp.name,
            '1', str(conf_count)
        ]
        subprocess.run(command)
        zrank_out_file = f'{tmp.name}.zr.out'
    
    # grab results from ranked file
    with open (zrank_out_file, 'r') as zrank_out:
        for line in zrank_out:
            pose, score = line.split('\t')
            pose = int(pose)
            score = float(score)

            pose_obj = ligase_obj.conformations[pose - 1]
            pose_obj.z_score = score

