import subprocess
from os import putenv
from pathlib import Path
from openbabel import pybel
from ..definitions import ROOT_DIR, CWD
from ..tools.logger import logger
from ..tools.structure_tools import pymol_combine

def rxdock_rescore(params):

    params['linker_scoring_folder'].mkdir(exist_ok=True)

    if params['minimize_protac']: dock_prm_file = 'minimize_score.prm'
    else: dock_prm_file = 'score.prm'

    # files
    folder = params['linker_scoring_folder']/f"protein_pose_{params['pose_number']}"
    combined_file = folder/'combined_protein.mol2'
    prep_protacs_file = folder/'protac_prep_confs.sdf'
    cavity_input_file = folder/'cavity.prm'
    scored_protacs_file = folder/'protac_scored_confs.sd'
    ##

    folder.mkdir(exist_ok=True)
    # set to rxdock that this will be a folder with definitions (minimize)
    putenv('RBT_HOME', str(Path(ROOT_DIR)/'inputs'))
    # make combined receptor+ligase file with pymol
    pymol_combine(params['receptor_obj_file'], params['file'], out_filename=combined_file)

    # prepare protac sdf file with openbabel
    sampled_confs = pybel.readfile('sdf', str(params['protac_pose']['file']))
    prep_confs = pybel.Outputfile('sdf', str(prep_protacs_file), overwrite=True)
    for mol in sampled_confs:
        mol.OBMol.DeleteNonPolarHydrogens()
        prep_confs.write(mol)

    # prepare cavity
    cavity_input = open(Path(ROOT_DIR)/'inputs'/'cavity.prm').read()
    cavity_input = cavity_input.replace('[title]', combined_file.stem)
    cavity_input = cavity_input.replace('[receptor]', combined_file.name)
    cavity_input = cavity_input.replace('[conformations]', prep_protacs_file.name)
    with open(cavity_input_file,'w+') as cavity_file:
        cavity_file.write(cavity_input)
    
    # run cavity generation
    subprocess.run(
        ['rbcavity', '-r', cavity_input_file.name, '-W'],
        cwd=Path(CWD)/folder, # we have to change working dirs for rxdock to work properly
        stdout=subprocess.DEVNULL
    )

    # run rxdock
    command = [
        'rbdock', '-i', prep_protacs_file.name, '-o', scored_protacs_file.stem,
        '-r', cavity_input_file.name,
        '-p', dock_prm_file
    ]
    subprocess.run(command, cwd=Path(CWD)/folder, stdout=subprocess.DEVNULL)

    # update dictionary with new scored file
    params['protac_pose']['scored_file'] = scored_protacs_file

    return(params)


def capture_rxdock_scores(params):

    from openbabel import pybel

    confs = pybel.readfile('sdf', str(params['protac_pose']['scored_file']))
    for conf in confs:
        conf_number = int(conf.data['Name'].split('_')[-1]) # name in sdf file is format "conf_X"
        score = float(conf.data['SCORE'])

        params['linker_confs'][conf_number]['rx_score'] = score

    return(params)