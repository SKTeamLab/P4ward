
if __name__ == '__main__':

    # prepare internals:
    from .config import config
    from .tools import run_tracker
    from .definitions import ROOT_DIR, CPT_FILE

    args = config.arg_parser(None)
    conf = config.make_config(args.config_file, ROOT_DIR)
    # load run tracker from previous run
    tracker = run_tracker.load_tracker(overwrite=conf.getboolean('general', 'overwrite'))

    # load pickle from previous run
    receptor, ligase = run_tracker.load_protein_objects(
        pickle_file=CPT_FILE,
        conf=conf,
        overwrite=conf.getboolean('general', 'overwrite')
    )


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~ start docking ~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # run megadock
    from .run import megadock

    megadock.prep_structures(
        receptor_obj=receptor,
        ligase_obj=ligase,
        choice=conf.getboolean('megadock', 'run_docking')
    )

    megadock.run_docking(
        program_path=conf.get('megadock', 'program_path'),
        receptor_file=receptor.mg_file,
        ligase_file=ligase.mg_file,
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

    # CHECKPOINT!
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase)

    megadock.zrank_rescore(
        ligase_obj=ligase,
        receptor_obj=receptor,
        zrank_path=conf.get('megadock', 'zrank_path'),
        run_docking_output_file=conf.get('megadock', 'run_docking_output_file'),
        choice=conf.getboolean('megadock', 'zrank_rescore'),
    )

    megadock.filter_poses(
        receptor_obj=receptor,
        ligase_obj=ligase,
        dist_cutoff=conf.getfloat('megadock','filter_dist_cutoff'),
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    # CHECKPOINT!
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase)

    megadock.cluster(
        pose_objects=ligase.active_confs(), # structures are now the saved confs attr
        clustering_cutoff=conf.getfloat('megadock','clustering_cutoff'),
        choice=conf.getboolean('megadock', 'cluster_poses')
    )


    #rank final protein poses
    from .analyse import rank

    rank.protein_poses(
        pose_objs=ligase.conformations,
        top_poses=conf.getint('protein_ranking', 'top_poses'),
        final_ranking_megadock_score=conf.getboolean('protein_ranking', 'final_ranking_megadock_score'),
        final_ranking_z_score=conf.getboolean('protein_ranking', 'final_ranking_z_score'),
        use_only_cluster_centroids=conf.getboolean('protein_ranking', 'use_only_cluster_centroids'),
        top_poses_from_centroids_only=conf.getboolean('protein_ranking', 'top_poses_from_centroids_only'),
    )

    rank.generate_protein_poses(
        pose_objs=ligase.conformations,
        poses=conf.get('protein_ranking', 'generate_poses'),
        generated_poses_folder=conf.get('protein_ranking', 'generated_poses_folder')
    )

    # CHECKPOINT!
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~ ligand sampling ~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    from .run import linker_sampling

    linker_sampling.rdkit_sampling(
        receptor_obj=receptor,
        ligase_obj=ligase,
        protac=conf.get('general', 'protac'),
        rdkit_number_of_confs=conf.getint('linker_sampling', 'rdkit_number_of_confs'),
        protac_poses_folder=conf.get('linker_sampling', 'protac_poses_folder'),
        rmsd_tolerance=conf.getfloat('linker_sampling', 'rdkit_pose_rmsd_tolerance'),
        choice=conf.getboolean('linker_sampling', 'rdkit_sampling'),
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~ end session ~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    from .analyse.make_summary import summary_csv

    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase)
    summary_csv(ligase.active_confs())