from ctypes import c_uint32, c_uint, c_ulong, c_char_p, c_double
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

def mod2pi(x):
    try:
        x = float(x)
    except:
        ValueError("Argument of mod2pi needs to be a float.")
    clibrebound.reb_tools_mod2pi.restype = c_double
    return clibrebound.reb_tools_mod2pi(c_double(x))

def M_to_f(e, M):
    try:
        e = float(e)
        M = float(M)
    except:
        ValueError("Arguments of M_to_f need to be floats.")
    clibrebound.reb_tools_M_to_f.restype = c_double
    return clibrebound.reb_tools_M_to_f(c_double(e), c_double(M))

def E_to_f(e, E):
    try:
        e = float(e)
        E = float(E)
    except:
        ValueError("Arguments of E_to_f need to be floats.")
    clibrebound.reb_tools_E_to_f.restype = c_double
    return clibrebound.reb_tools_E_to_f(c_double(e), c_double(E))

def M_to_E(e, M):
    try:
        e = float(e)
        M = float(M)
    except:
        ValueError("Arguments of M_to_E need to be floats.")
    clibrebound.reb_tools_M_to_E.restype = c_double
    return clibrebound.reb_tools_M_to_E(c_double(e), c_double(M))
