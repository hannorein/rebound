from ctypes import c_uint32, c_char_p
from . import clibrebound

def hash(string):
    clibrebound.reb_tools_hash.restype = c_uint32
    return clibrebound.reb_tools_hash(c_char_p(string.encode('utf-8')))

