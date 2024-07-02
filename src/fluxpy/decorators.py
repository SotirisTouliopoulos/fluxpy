import functools

def print_doc_on_typeerror(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except TypeError as e:
            print(f"TypeError occurred in function '{func.__name__}': {e}")
            if func.__doc__:
                print("Function description:")
                print(func.__doc__)
            raise
    return wrapper