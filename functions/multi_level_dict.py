"""
Description:
    function to create multi-level dictionaries used extensively throughout.
    
Author
    Ravi kumar verma
"""

from numbers import Number


class mutileveldict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

    def __add__(self, x):
        """ override addition for numeric types when self is empty """
        if not self and isinstance(x, Number):
            return x
        raise ValueError

    def __sub__(self, x):
        if not self and isinstance(x, Number):
            return -1 * x
        raise ValueError
