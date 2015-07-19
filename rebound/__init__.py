# -*- coding: utf-8 -*-

"""
REBOUND Module. 

"""
#Make changes for python 2 and 3 compatibility
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x

from ctypes import *
import ctypes.util
import os
_pymodulespath = os.path.dirname(__file__)
try:
    clibrebound = CDLL(_pymodulespath+"/../librebound.so", RTLD_GLOBAL)
except:
    print("Cannot find library 'librebound.so'.")
    raise

def build_str():
    return str(c_char_p.in_dll(clibrebound, "reb_build_str").value)

from .librebound import Simulation
   
def create_simulation():
    clibrebound.reb_create_simulation.restype = Simulation
    return clibrebound.reb_create_simulation()
    
    

from .particle import Particle
from .particle import Orbit
from .interruptible_pool import InterruptiblePool

__all__ = ["Simulation", "Orbit", "Particle"]
