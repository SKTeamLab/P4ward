import subprocess
import os
from ..tools import decorators
from ..tools.script_tools import logit
from ..definitions import CWD


@decorators.user_choice
def run(program_path, receptor_file, ligase_file, num_threads,
        num_predictions, num_predictions_per_rotation,
        logger=None):
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

    logit.info('Running megadock...')
    with open('megadock_run.out', 'w+') as outfile: 
        subprocess.run(command, stdout=outfile)
    logit.info('Done.')


@decorators.user_choice
def filter_poses(receptor, receptor_ligand_resnum, ligase, ligase_ligand_resnum,
                 receptor_ligand_chain, ligase_ligand_chain,
                 dist_cutoff,
                 #output_file, TODO add output support
                 #output_filtered_file,
                 logger=None):
    # TODO ver se aqui vai passar a receber o obj pronto do biopython
    """
    Use Biopython to filter the megadock poses which satisfy a
    distance cutoff for both binding sites
    """
    from ..tools.structure_tools import load_biopython_structures
    from ..tools.structure_tools import structure_proximity
    from ..tools.classes import ProteinPose

    logit.info(f'Filtering megadock poses with cuttoff {dist_cutoff}')

    # first get receptor, and its ligand, they will not change
    _, receptor_ligand_struct = load_biopython_structures(
        protein=receptor,
        protein_ligand_chain=receptor_ligand_chain,
        protein_ligand_resnum=receptor_ligand_resnum,
        logger=logger
    )

    output_file = 'megadock.out'
    output_filtered_file = 'megadock-filtered.out'
    output = open(output_file, 'r')
    output_filtered = open(output_filtered_file, 'a+')

    # make a folder for the docked structures
    output_path = os.path.join(CWD, 'protein_docking')
    if os.path.isdir(output_path):
        mssg = f'Docked structure folder {output_path} exists.'
        try: logit.info(mssg)
        except: pass
    else:
        os.mkdir(output_path)

    count = 0
    for line in output:

        if len(line.split('\t')) == 7:
            count += 1

            decoy_name = os.path.join(output_path, f'decoy{count}.pdb')
            decoygen_command = [
                'decoygen', decoy_name,
                ligase,
                output_file, str(count)
            ]
            subprocess.run(decoygen_command, stdout=subprocess.DEVNULL)

            _, ligase_ligand_struct = load_biopython_structures(
                protein=decoy_name,
                protein_ligand_chain=ligase_ligand_chain,
                protein_ligand_resnum=ligase_ligand_resnum,
                logger=logger
            )

            distance, proximity = structure_proximity(
                receptor_ligand_struct,
                ligase_ligand_struct,
                dist_cutoff=dist_cutoff
            )

            structures = []
            if proximity:
                output_filtered.write(line)
                structures.append(decoy_name)
                logit.info(f'Saving filtered pose {count}, with distance of {distance} between ligands.')
            else:
                logit.info(f'Ignoring pose {count}, with distance of {distance} between ligands.')
                os.remove(decoy_name)

        else: 
            output_filtered.write(line)
    
    logit.info(f'Saved filtered megadock poses to file {output_filtered_file}')

    output.close()
    output_filtered.close()

    return(structures, output_filtered_file)


@decorators.user_choice
def cluster(structure_list, clustering_cutoff, logger=None):
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



