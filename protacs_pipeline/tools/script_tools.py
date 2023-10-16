import os
from shutil import copy
from pathlib import Path
from ..definitions import CWD, ROOT_DIR


# def create_folder(folder_path):
#     """
#     creates a folder if it does not exist, does nothing if it does
#     logs what happened
#     """

#     if os.path.isdir(folder_path):
#         logger.info(f'Folder {folder_path} exists.')
#     else:
#         os.mkdir(folder_path)
#         logger.info(f'Created {folder_path}')


def write_default_conf():
    """
    replicate default.ini into working dir for user reference.
    """

    src = Path(ROOT_DIR)/'config'/'default.ini'
    copy(src=str(src), dst=CWD)

    print("Copied default.ini configuration file into current working directory.")