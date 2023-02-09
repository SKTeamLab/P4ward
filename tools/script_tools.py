import logging
from ..config import config
import os
import datetime as dt

def set_logging(PIPELINE_NAME):
    """
    Configure the logger for the whole pipeline
    """

    now = dt.datetime.now()

    logfile_name = os.path.join(os.getcwd(), PIPELINE_NAME + '.log')
    with open(logfile_name, 'w+') as logfile:
        logfile.write(
            f"""
-----------------------------------
PROTACS PIPELINE
Log created at {now.strftime("%Y-%m-%d, %H:%M:%S")}
-----------------------------------
\n"""
        )

    file_format = logging.Formatter(
        '%(asctime)s > %(levelname)s - %(message)s',
        "%H:%M:%S"
    )

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(file_format)
    file_handler = logging.FileHandler(logfile_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_format)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return(logger)


class logit:
    """
    Log message to logger and, if logger does not exist (when __name__ != '__main__'),
    just print message to the terminal
    """

    def __init__(self) -> None:
        pass

    def critical(log_message):
        try:
            logger.critical(log_message)
        except:
            print(f'CRITICAL > {log_message}')
            
    def error(log_message):
        try:
            logger.error(log_message)
        except:
            print(f'ERROR > {log_message}')

    def warning(log_message):
        try:
            logger.warning(log_message)
        except:
            print(f'WARNING > {log_message}')

    def info(log_message):
        try:
            logger.info(log_message)
        except:
            print(f'INFO > {log_message}')

    def debug(log_message):
        try:
            logger.debug(log_message)
        except:
            print(f'DEBUG > {log_message}')

    def notset(log_message):
        try:
            logger.notset(log_message)
        except:
            print(f'NOTSET > {log_message}')

  

