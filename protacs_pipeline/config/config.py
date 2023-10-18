import configparser
import argparse
from shutil import copy
from pathlib import Path
from ..definitions import CWD, ROOT_DIR


def arg_parser(arguments):
    """
    Make argument parser
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config_file",
        type=Path,
        help='Path to the .ini file where the user configuration is stored'
    )
    parser.add_argument(
        "--write_default",
        action='store_true',
        help='Copy full configuration file and quit the program.'
    )
    parser.add_argument(
        "--benchmark",
        action='store_true',
        help='Benchmark a modelling run.'
    )
    parser.add_argument(
        "--ref_ligase",
        type=Path,
        help="Path to the reference pose for the E3 ligase."
    )

    if arguments == None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(arguments)
    return(args)


def make_config(config_file, ROOT_DIR):
    """
    read default config + user config, and return object
    """

    conf = configparser.ConfigParser()
    conf.read(ROOT_DIR/'config'/'default.ini')
    conf.read(config_file)

    return(conf)


def write_default_conf():
    """
    replicate default.ini into working dir for user reference.
    """

    src = Path(ROOT_DIR)/'config'/'default.ini'
    copy(src=str(src), dst=CWD)

    print("Copied default.ini configuration file into current working directory.")