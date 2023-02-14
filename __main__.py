
if __name__ == '__main__':

    # prepare internals:
    from .config import config
    from .tools.logger import logger
    from .tools.script_tools import load_protein_objects, save_protein_objects
    from .definitions import ROOT_DIR, PICKLE_FILE

    args = config.arg_parser(None)
    conf = config.make_config(args.config_file, ROOT_DIR)

    # load pickle from previous run
    receptor, ligase = load_protein_objects(
        pickle_file=PICKLE_FILE,
        conf=conf,
        overwrite=conf.getboolean('general', 'overwrite')
    )
    # load the biopython objects as attributes for receptor and ligand:
    receptor.protein_struct = receptor.get_protein_struct()
    receptor.ligand_struct = receptor.get_ligand_struct()
    ligase.protein_struct = ligase.get_protein_struct()
    ligase.ligand_struct = ligase.get_ligand_struct()


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~ start running ~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # run megadock
    from .run import megadock

    megadock.run(
        program_path=conf.get('megadock', 'program_path'),
        receptor_file=receptor.file,
        ligase_file=ligase.file,
        num_threads=conf.get('megadock', 'num_threads'),
        run_docking_output_file=conf.get('megadock', 'run_docking_output_file'),
        num_predictions=conf.get('megadock', 'num_predictions'),
        num_predictions_per_rotation=conf.get('megadock', 'num_predictions_per_rotation'),
        log_file=conf.get('megadock', 'run_docking_log_file'),
        choice=conf.getboolean('megadock', 'run_docking')
    )

    megadock.capture_scores(
        ligase_obj=ligase,
        run_docking_output_file=conf.get('megadock', 'run_docking_output_file'),
        choice=conf.getboolean('megadock', 'run_docking') # runs if docking is run
    )

    megadock.generate_poses(
        ligase_obj=ligase,
        run_docking_output_file=conf.get('megadock', 'run_docking_output_file'),
        docked_poses_folder=conf.get('megadock', 'docked_poses_folder'),
        choice=conf.getboolean('megadock', 'generate_all_poses'),
    )

    megadock.filter_poses(
        receptor_obj=receptor,
        ligase_obj=ligase,
        dist_cutoff=conf.getfloat('megadock','filter_dist_cutoff'),
        output_file=conf.get('megadock', 'run_docking_output_file'),
        output_filtered_file=conf.get('megadock', 'filter_poses_output_file'),
        generate_all_poses=conf.getboolean('megadock', 'generate_all_poses'),
        docked_poses_folder=conf.get('megadock', 'docked_poses_folder'),
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    megadock.cluster(
        pose_objects=ligase.active_confs(), # structures are now the saved confs attr
        clustering_cutoff=conf.getfloat('megadock','clustering_cutoff'),
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~ end session ~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    save_protein_objects(
        receptor_obj=receptor,
        ligase_obj=ligase,
        pickle_file=PICKLE_FILE
    )