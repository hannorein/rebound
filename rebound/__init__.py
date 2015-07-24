# -*- coding: utf-8 -*-

"""
REBOUND Module. 

"""
# Make changes for python 2 and 3 compatibility
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x

# Find suffix
import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Import shared library
from ctypes import *
clibrebound = cdll.LoadLibrary("librebound"+suffix)

def build_str():
    return str(c_char_p.in_dll(clibrebound, "reb_build_str").value)


# Exceptions    
class SimulationError(Exception):  
    pass

class Encounter(Exception):
    pass

class Escape(Exception):
    pass

class NoParticles(Exception):
    pass
    


from .simulation import Simulation
from .particle import Particle
from .particle import Orbit
from .interruptible_pool import InterruptiblePool

__all__ = ["Simulation", "Orbit", "Particle", "SimulationError", "Encounter", "Escape", "NoParticles"]
