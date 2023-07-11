from pathlib import Path
from os import putenv
import subprocess
from ..definitions import ROOT_DIR, CWD
from ..tools import decorators
from ..tools.logger import logger
from ..tools.script_tools import create_folder
from ..tools.structure_tools import pymol_combine


@decorators.user_choice
@decorators.track_run
def rxdock_rescore(
                        receptor_obj,
                        ligase_obj,
                        minimize,
):
    
    from openbabel import pybel

    # we will operate on all the pose objs that have at least an active linker conf
    pose_objs = [i for i in ligase_obj.active_confs() if i.protac_pose.active and len(i.protac_pose.linker_confs) > 0]
    
    create_folder('./ligand_scoring')

    if minimize: dock_prm_file = 'minimize_score.prm'
    else: dock_prm_file = 'score.prm'

    for pose_obj in pose_objs:

        # files
        folder = Path(f'./ligand_scoring/protein_pose_{pose_obj.pose_number}')
        combined_file = folder/'combined_protein.mol2'
        prep_protacs_file = folder/'protac_prep_confs.sdf'
        cavity_input_file = folder/'cavity.prm'
        scored_protacs_file = folder/'protac_scored_confs.sd'
        ##

        create_folder(folder)
        # set to rxdock that this will be a folder with definitions (minimize)
        putenv('RBT_HOME', str(Path(ROOT_DIR)/'inputs'))
        # make combined receptor+ligase file with pymol
        pymol_combine(receptor_obj.file, pose_obj.file, out_filename=combined_file)

        # prepare protac sdf file with openbabel
        sampled_confs = pybel.readfile('sdf', pose_obj.protac_pose.file)
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
            cwd=Path(CWD)/folder # we have to change working dirs for rxdock to work properly
        )
        logger.info("Ran cavity detection for rxdock rescoring.")

        # run rxdock
        command = [
            'rbdock', '-i', prep_protacs_file.name, '-o', scored_protacs_file.stem,
            '-r', cavity_input_file.name,
            '-p', dock_prm_file
        ]
        subprocess.run(command, cwd=Path(CWD)/folder)
        logger.info("Ran rxdock rescoring.")

        # update object
        pose_obj.protac_pose.scored_file = scored_protacs_file


@decorators.user_choice
# @decorators.track_run
def capture_rxdock_scores(ligase_obj):

    from openbabel import pybel

    # get the protac poses that have a scored file by rxdock
    protac_poses = [i.protac_pose for i in ligase_obj.active_confs() if hasattr(i.protac_pose, 'scored_file')]

    for protac_pose in protac_poses:

        confs = pybel.readfile('sdf', str(protac_pose.scored_file))
        for conf in confs:
            conf_number = conf.data['Name'].split('_')[-1] # name in sdf file is format "conf_X"
            score = conf.data['SCORE']

            linker_conf = [i for i in protac_pose.linker_confs if i.conf_number == int(conf_number)][0]
            linker_conf.rx_score = score
