


if __name__ == '__main__':

    from sys import exit
    from pathlib import Path
    from .config import config
    from .tools import run_tracker
    from .tools.script_tools import write_default_conf
    from .definitions import ROOT_DIR, CPT_FILE

    # prepare internals:
    args = config.arg_parser(None)

    if args.write_default:
        write_default_conf()
        exit(0)

    conf = config.make_config(args.config_file, ROOT_DIR)
    # load run tracker from previous run
    tracker = run_tracker.load_tracker(overwrite=conf.getboolean('general', 'overwrite'))

    # load pickle from previous run
    receptor, ligase, protac = run_tracker.load_run_objects(
        pickle_file=CPT_FILE,
        conf=conf,
        overwrite=conf.getboolean('general', 'overwrite')
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~ prep proteins ~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
   
    from .run import md

    md.fix_proteins(
        receptor,
        ligase,
        fixed_suffix='_fixed',
        ignore_extremities=conf.get('protein_prep', 'pdbfixer_ignore_extremities'),
        ph=conf.getfloat('protein_prep', 'pdbfixer_ph'),
        choice=conf.getboolean('protein_prep', 'pdbfixer')
    )

    md.minimize_proteins(
        receptor,
        ligase,
        minimized_suffix='_minim',
        maxiterations=conf.getint('protein_prep', 'minimize_maxiter'),
        choice=conf.getboolean('protein_prep', 'minimize')
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~ start docking ~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    # run megadock
    from .run import megadock
    from .tools import structure_tools

    megadock.prep_structures(
        receptor_obj=receptor,
        ligase_obj=ligase,
        choice=conf.getboolean('megadock', 'run_docking')
    )

    megadock.run_docking(
        program_path=conf.get('program_paths', 'megadock'),
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
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase, protac_obj=protac)

    megadock.zrank_rescore(
        ligase_obj=ligase,
        receptor_obj=receptor,
        zrank_path=conf.get('megadock', 'zrank_path'),
        run_docking_output_file=conf.get('megadock', 'run_docking_output_file'),
        choice=conf.getboolean('megadock', 'zrank_rescore'),
    )

    structure_tools.get_protac_dist_cuttoff(
        protac_obj=protac,
        reclig_file=receptor.lig_file,
        liglig_file=ligase.lig_file,
        dist_cutoff=conf.get('megadock', 'filter_dist_cutoff'),
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    megadock.filter_poses(
        receptor_obj=receptor,
        ligase_obj=ligase,
        protac_obj=protac,
        dist_cutoff=protac.dist_cutoff,
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    # CHECKPOINT!
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase, protac_obj=protac)

    megadock.cluster(
        pose_objects=ligase.active_confs(), # structures are now the saved confs attr
        clustering_cutoff=conf.getfloat('megadock','clustering_cutoff'),
        choice=conf.getboolean('megadock', 'cluster_poses')
    )


    #rank final protein poses
    from .analyse import rank

    rank.protein_poses(
        ligase_obj=ligase,
        top_poses=conf.getint('protein_ranking', 'top_poses'),
        final_ranking_megadock_score=conf.getboolean('protein_ranking', 'final_ranking_megadock_score'),
        final_ranking_z_score=conf.getboolean('protein_ranking', 'final_ranking_z_score'),
        cluster_rep=conf.get('protein_ranking', 'cluster_rep'),
        rank_cluster_reps_only=conf.getboolean('protein_ranking', 'rank_cluster_reps_only'),
    )

    # CHECKPOINT!
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase, protac_obj=protac)


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~ ligand sampling ~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    from .run import protac_sampling

    # generate pdbs
    rank.generate_protein_poses(
        pose_objs=ligase.conformations,
        poses=conf.get('protein_ranking', 'generate_poses'),
        altlocA = conf.getboolean('protein_ranking', 'generate_poses_altlocA'),
        generated_poses_folder=Path(conf.get('protein_ranking', 'generated_poses_folder'))
    )


    # CHECKPOINT!
    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase, protac_obj=protac)

    # run protac sampling

    protac_sampling.protac_sampling(
        receptor_obj=receptor,
        ligase_obj=ligase,
        protac_obj=protac,
        extend_flexible_small_linker=conf.getboolean('linker_sampling', 'extend_flexible_small_linker'),
        neighbour_number=conf.getint('linker_sampling', 'extend_neighbour_number'),
        min_linker_length=conf.getint('linker_sampling', 'min_linker_length'),
        rdkit_number_of_confs=conf.getint('linker_sampling', 'rdkit_number_of_confs'),
        protac_poses_folder=Path(conf.get('linker_sampling', 'protac_poses_folder')),
        rmsd_tolerance=conf.getfloat('linker_sampling', 'rdkit_pose_rmsd_tolerance'),
        time_tolerance=conf.getint('linker_sampling', 'rdkit_time_tolerance'),
        extend_top_poses_sampled=conf.getboolean('linker_sampling', 'extend_top_poses_sampled'),
        extend_top_poses_score=conf.getboolean('linker_sampling', 'extend_top_poses_score'),
        extend_top_poses_energy=conf.getboolean('linker_sampling', 'extend_top_poses_energy'),
        linker_scoring_folder=Path(conf.get('linker_ranking','linker_scoring_folder')),
        minimize_protac=conf.getboolean('linker_ranking','rxdock_minimize'),
        num_parallel_procs=conf.getint('general', 'num_processors'),
        choice=conf.getboolean('linker_sampling', 'rdkit_sampling')
    )

    # rank all sampled protacs
    # rank.protac_conformations(
    #     protac_poses=protac.active_poses()
    # )


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~ end session ~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    from .analyse import summaries

    run_tracker.save_protein_objects(receptor_obj=receptor, ligase_obj=ligase, protac_obj=protac)
    summaries.summary_csv([i for i in ligase.conformations if i.top])
    summaries.chimerax_view(
        receptor_obj=receptor,
        pose_objs=[i for i in ligase.conformations if i.top],
        generated_poses_folder=conf.get('protein_ranking', 'generated_poses_folder'),
        protac_poses_folder=conf.get('linker_sampling', 'protac_poses_folder')
    )