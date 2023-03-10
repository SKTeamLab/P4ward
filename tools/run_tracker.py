import pickle
import os
from .logger import logger
from . import classes
from ..definitions import TRACKER_FILE, CPT_FILE


def load_tracker(overwrite):
    """
    create/retrieve the dictionary that tracks which function was run
    and dump it to the pickle object. Should be run only once in the
    beginning of the run.
    """
    if not overwrite and os.path.isfile(TRACKER_FILE):

        logger.info('Retrieving previous run steps.')
        with open(TRACKER_FILE, 'rb') as pic:
            run_tracker = pickle.load(pic)
    else:
        logger.info('Not retrieving previous run steps.')
        run_tracker = []
    
    with open(TRACKER_FILE, 'wb+') as pic:
        pickle.dump(run_tracker, pic)


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
            lig_file=conf.get('general','receptor_ligand')
        )
        ligase_obj = classes.Protein(
            ptn_type='ligase',
            file=conf.get('general', 'ligase'),
            lig_file=conf.get('general','ligase_ligand')
        )

    return(receptor_obj, ligase_obj)


def save_protein_objects(receptor_obj, ligase_obj):
    """
    Save receptor and ligase objects into pickle file
    """

    with open(CPT_FILE, 'wb+') as pic:
        pickle.dump((receptor_obj, ligase_obj), pic)
    
    logger.info('Saved information from this run.')
