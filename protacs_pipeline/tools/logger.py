import logging
import datetime as dt
from ..definitions import PIPELINE_NAME, CWD


def set_logging():
    """
    Configure the logger for the whole pipeline
    """

    now = dt.datetime.now()

    logfile_name = CWD/(PIPELINE_NAME + '.log')
    with open(logfile_name, 'w+') as logfile:
        logfile.write(
            f"""
-----------------------------------
{PIPELINE_NAME}
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
  

logger = set_logging()

  

