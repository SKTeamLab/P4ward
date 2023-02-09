import configparser
import argparse
import os

# def get_paths():
#     """
#     get root dir and cwd global variables
#     """

#     ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
#     CWD = os.getcwd()

#     return(ROOT_DIR, CWD)


def arg_parser(arguments):
    """
    Make argument parser
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config_file",
        type=str,
        help='Path to the .ini file where the user configuration is stored'
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
    conf.read(os.path.join(ROOT_DIR, 'config', 'default.ini'))
    conf.read(config_file)

    return(conf)
