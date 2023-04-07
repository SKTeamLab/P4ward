import functools
import pickle
from ..tools.logger import logger
from ..definitions import TRACKER_FILE, CPT_FILE

def user_choice(func):
    """
    For each (applicable) function, there will be choice made by
    the user of whether to run it or not. This decorator makes it
    easy to implement this for every single one of those functions.
    """

    @functools.wraps(func)
    def wrapper(*args, choice=True, **kwargs):
        if choice:
            func(*args, **kwargs)
        else:
            pass
    
    return(wrapper)


def track_run(func):
    """
    Track if a function has been called by checking if its name is in run_tracker
    if not, run it, then add its name to the list
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        
        with open(TRACKER_FILE, 'rb') as pic:
            tracker = pickle.load(pic)
        
        func_path = f"{func.__module__}.{func.__name__}"

        if func_path not in tracker:

            func(*args, **kwargs)

            tracker.append(func_path)
            with open(TRACKER_FILE, 'wb+') as pic:
                pickle.dump(tracker, pic)
            
        else:
            logger.debug(f'Skipping run for function {func_path}')
    return(wrapper)
