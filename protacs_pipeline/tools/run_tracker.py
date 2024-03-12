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
    if not overwrite and TRACKER_FILE.is_file():

        logger.info('Retrieving previous run steps.')
        with open(TRACKER_FILE, 'rb') as pic:
            run_tracker = pickle.load(pic)
    else:
        logger.info('Not retrieving previous run steps.')
        run_tracker = []
    
    with open(TRACKER_FILE, 'wb+') as pic:
        pickle.dump(run_tracker, pic)


def load_run_objects(pickle_file, conf, overwrite=False):
    """
    load all objects either from previously saved pickle or
    from scratch when first time running or overwrite chosen by user
    """

    def parse_protacs(protacs_file):

        protac_objs = []
        with open(protacs_file, 'r') as smiles_file:
            count = 0
            for line in smiles_file:
                count += 1
                line = line.strip().split()
                if len(line) == 2:
                    smiles = line[0]
                    name = line[1]
                elif len(line) == 1:
                    smiles = line[0]
                    name = f"protac{count}"
                protac_obj = classes.Protac(smiles=smiles, name=name, number=count)
                protac_objs.append(protac_obj)
        return(protac_objs)
                
                
    if not overwrite and pickle_file.is_file():
        logger.info('Retrieving information from previous run.')
        with open(pickle_file, 'rb') as pic:
            receptor_obj, ligase_obj, protac_objs = pickle.load(pic)
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
        protac_objs = parse_protacs(conf.get('general', 'protacs'))


    return(receptor_obj, ligase_obj, protac_objs)


def save_protein_objects(receptor_obj, ligase_obj, protac_objs):
    """
    Save all objects into pickle file
    """

    with open(CPT_FILE, 'wb+') as pic:
        pickle.dump((receptor_obj, ligase_obj, protac_objs), pic)
    
    logger.info('Saved information from this run.')
