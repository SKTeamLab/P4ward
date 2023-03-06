import os
from .logger import logger

def create_folder(folder_path):
    """
    creates a folder if it does not exist, does nothing if it does
    logs what happened
    """

    if os.path.isdir(folder_path):
        logger.info(f'Folder {folder_path} exists.')
    else:
        os.mkdir(folder_path)
        logger.info(f'Created {folder_path}')