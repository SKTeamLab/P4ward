import configparser
import argparse
from pathlib import Path


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
