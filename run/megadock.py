import subprocess
import os
from ..tools import decorators
from ..tools.logger import logger
from ..definitions import CWD

@decorators.user_choice
def run(program_path, receptor_file, ligase_file, num_threads,
        num_predictions, num_predictions_per_rotation):
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
        '-o', 'megadock.out'
    ]

    logger.info('Running megadock...')
    with open('megadock_run.out', 'w+') as outfile: 
        subprocess.run(command, stdout=outfile)
    logger.info('Done.')


@decorators.user_choice
def filter_poses(receptor_obj, ligase_obj,
                 dist_cutoff,
                 #output_file, TODO add output support
                 #output_filtered_file,
                 ):
    """
    Use Biopython to filter the megadock poses which satisfy a
    distance cutoff for both binding sites. Takes the a Protein object
    and handles the rest by accessing its attributes
    """

    from ..tools.structure_tools import structure_proximity
    from ..tools.classes import ProteinPose

    logger.info(f'Filtering megadock poses with cuttoff {dist_cutoff}')

    output_file = 'megadock.out'
    output_filtered_file = 'megadock-filtered.out'
    output = open(output_file, 'r')
    output_filtered = open(output_filtered_file, 'a+')

    # make a folder for the docked structures
    output_path = os.path.join(CWD, 'protein_docking')
    if os.path.isdir(output_path):
        mssg = f'Docked structure folder {output_path} exists.'
        logger.info(mssg)
    else:
        os.mkdir(output_path)

    count = 0
    for line in output:

        if len(line.split('\t')) == 7:
            count += 1

            decoy_name = os.path.join(output_path, f'decoy{count}.pdb')
            decoygen_command = [
                'decoygen', decoy_name,
                ligase_obj.file,
                output_file, str(count)
            ]
            subprocess.run(decoygen_command, stdout=subprocess.DEVNULL)

            ligase_pose_obj = ProteinPose(parent=ligase_obj, file=decoy_name)

            distance, proximity = structure_proximity(
                receptor_obj.ligand_struct,
                ligase_pose_obj.get_ligand_struct(),
                dist_cutoff=dist_cutoff
            )

            if proximity:
                output_filtered.write(line)
                ligase_pose_obj.active = True
                logger.info(f'Saving filtered pose {count}, with distance of {distance} between ligands.')
            else:
                ligase_pose_obj.active = False
                ligase_pose_obj.file = None
                logger.info(f'Ignoring pose {count}, with distance of {distance} between ligands.')
                os.remove(decoy_name)

        else: 
            output_filtered.write(line)
    
    logger.info(f'Saved filtered megadock poses to file {output_filtered_file}')

    output.close()
    output_filtered.close()

    return(output_filtered_file)


@decorators.user_choice
def cluster(structure_list, clustering_cutoff):
    """
    cluster megadock docked poses using a user-specified cutoff
    """

    from ..tools.structure_tools import load_biopython_structures
    from sklearn.cluster import AgglomerativeClustering
    from ..tools.structure_tools import get_rmsd
    import numpy as np
    
    # get first ligase as reference:
    reference_structure_struct = load_biopython_structures(protein=structure_list[0])

    rmsd_values = {'name':[], 'rmsd':[]}
    # load receptor object
    for ligase in structure_list[10:]:
        ligase_name = os.path.basename(ligase)
        ligase_struct = load_biopython_structures(protein=ligase, logger=logger)

        rmsd = get_rmsd(reference_structure_struct, ligase_struct, ca=True)
        rmsd_values['name'].append(ligase_name)
        rmsd_values['rmsd'].append(rmsd)
    
    rmsd_points = np.asarray(rmsd_values['rmsd']).reshape(-1, 1)

    clustering = AgglomerativeClustering(
        n_clusters=None,
        linkage='average',
        distance_threshold=clustering_cutoff
    ).fit(rmsd_points)

    return(list(clustering.labels_))



