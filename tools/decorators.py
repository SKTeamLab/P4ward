import functools


def user_choice(func):
    """
    For each (applicable) function, there will be choice made by
    the user of whether to run it or not. This decorator makes it
    easy to implement this for every single one of those functions.
    """

    @functools.wraps(func)
    def wrapper(*args, choice=False, **kwargs):
        if choice:
            return func(*args, **kwargs)
        else:
            pass
    
    return(wrapper)
