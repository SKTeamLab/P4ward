import pickle
import os
from .logger import logger
from . import classes

def load_protein_objects(pickle_file, conf, overwrite=False):
    """
    load receptor and ligase objects either from previously saved pickle or
    from scratch when first time running or overwrite chosen by user
    """
    
    if not overwrite and os.path.isfile(pickle_file):
        logger.info('Retrieving information from previous run.')
        with open(pickle_file, 'rb') as pic:
            receptor_obj, ligase_obj = pickle.load(pic)
    else:
        logger.info('No previous data retrieved.')
        receptor_obj = classes.Protein(
            ptn_type='receptor',
            file=conf.get('general', 'receptor'),
            lig_chain=conf.get('general','receptor_ligand_chain'),
            lig_resnum=conf.getint('general','receptor_ligand_resnum')
        )
        ligase_obj = classes.Protein(
            ptn_type='ligase',
            file=conf.get('general', 'ligase'),
            lig_chain=conf.get('general','ligase_ligand_chain'),
            lig_resnum=conf.getint('general','ligase_ligand_resnum')
        )

    return(receptor_obj, ligase_obj)


def save_protein_objects(receptor_obj, ligase_obj, pickle_file):
    """
    Save receptor and ligase objects into pickle file
    """

    with open(pickle_file, 'wb+') as pic:
        pickle.dump((receptor_obj, ligase_obj), pic)
    
    logger.info('Saved information from this run.')
