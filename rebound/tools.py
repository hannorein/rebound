from ctypes import c_uint32, c_uint, c_ulong, c_char_p
from . import clibrebound
import sys

def hash(key):
    hash_types = c_uint32, c_uint, c_ulong
    PY3 = sys.version_info[0] == 3
    if PY3:
        string_types = str,
        int_types = int,
    else:
        string_types = basestring,
        int_types = int, long,
   
    if isinstance(key, int_types):
        return c_uint32(key)
    elif isinstance(key, hash_types):
        return key
    elif isinstance(key, string_types):
        clibrebound.reb_hash.restype = c_uint32
        return c_uint32(clibrebound.reb_hash(c_char_p(key.encode('ascii'))))
    else:
        raise AttributeError("Need to pash hash an integer or string.")

