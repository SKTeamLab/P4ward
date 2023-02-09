if __name__ == '__main__':

    # prepare internals:
    from .config import config
    from .tools import script_tools
    from .definitions import ROOT_DIR, CWD

    args = config.arg_parser(None)
    conf = config.make_config(args.config_file, ROOT_DIR)
    logger = script_tools.set_logging('protacs_pipeline')

    # run megadock
    from .run import megadock

    ## for now, since no other steps prior:
    receptor = conf.get('general', 'receptor')
    ligase = conf.get('general', 'ligase')

    megadock.run(
        program_path=conf.get('megadock', 'program_path'),
        receptor=receptor,
        ligase=ligase,
        num_threads=conf.get('megadock', 'num_threads'),
        num_predictions=conf.get('megadock', 'num_predictions'),
        num_predictions_per_rotation=conf.get('megadock', 'num_predictions_per_rotation'),
        logger=logger,
        choice=conf.getboolean('megadock', 'run_docking')
    )

    structures, output = megadock.filter_poses(
        receptor=receptor,
        ligase=ligase,
        receptor_ligand_resnum=conf.getint('general','receptor_ligand_resnum'),
        ligase_ligand_resnum=conf.getint('general','ligase_ligand_resnum'),
        receptor_ligand_chain=conf.get('general','receptor_ligand_chain'),
        ligase_ligand_chain=conf.get('general','ligase_ligand_chain'),
        dist_cutoff=conf.getfloat('megadock','filter_dist_cutoff'),
        logger=logger,
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    megadock.cluster(
        structure_list=structures, # structures are now `structures` from megadock.filter_poses
        clustering_cutoff=conf.getfloat('megadock','clustering_cutoff'),
        logger=logger,
        choice=conf.getboolean('megadock', 'filter_poses')
    )

