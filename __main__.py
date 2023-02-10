if __name__ == '__main__':

    # prepare internals:
    from .config import config
    from .tools.logger import logger
    from .tools import classes
    from .definitions import ROOT_DIR, CWD

    args = config.arg_parser(None)
    conf = config.make_config(args.config_file, ROOT_DIR)

    # initialise receptor and ligase objects with the definitions
    # from the configuration file:
    receptor = classes.Protein(
        type='receptor',
        file=conf.get('general', 'receptor'),
        lig_chain=conf.get('general','receptor_ligand_chain'),
        lig_resnum=conf.getint('general','receptor_ligand_resnum')
    )
    ligase = classes.Protein(
        type='ligase',
        file=conf.get('general', 'ligase'),
        lig_chain=conf.get('general','ligase_ligand_chain'),
        lig_resnum=conf.getint('general','ligase_ligand_resnum')
    )
    # we can load their biopython objects as attributes for receptor and ligand:
    receptor.protein_struct = receptor.get_protein_struct()
    receptor.ligand_struct = receptor.get_ligand_struct()
    ligase.protein_struct = ligase.get_protein_struct()
    ligase.ligand_struct = ligase.get_ligand_struct()

    # run megadock
    from .run import megadock

    megadock.run(
        program_path=conf.get('megadock', 'program_path'),
        receptor_file=receptor.file,
        ligase_file=ligase.file,
        num_threads=conf.get('megadock', 'num_threads'),
        num_predictions=conf.get('megadock', 'num_predictions'),
        num_predictions_per_rotation=conf.get('megadock', 'num_predictions_per_rotation'),
        choice=conf.getboolean('megadock', 'run_docking')
    )

    output = megadock.filter_poses(
        receptor_obj=receptor,
        ligase_obj=ligase,
        dist_cutoff=conf.getfloat('megadock','filter_dist_cutoff'),
        choice=conf.getboolean('megadock', 'filter_poses')
    )

    # megadock.cluster(
    #     structure_list=structures, # structures are now `structures` from megadock.filter_poses
    #     clustering_cutoff=conf.getfloat('megadock','clustering_cutoff'),
    #     choice=conf.getboolean('megadock', 'filter_poses')
    # )

