# -*- coding: utf-8 -*-
"""An N-body integrator package for python."""
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
import os
import warnings
pymodulepath = os.path.dirname(__file__)
from ctypes import cdll, c_char_p
__libpath__ = pymodulepath+"/../librebound"+suffix
clibrebound = cdll.LoadLibrary(__libpath__)

# Version
__version__ = c_char_p.in_dll(clibrebound, "reb_version_str").value.decode('ascii')

# Build
__build__ = c_char_p.in_dll(clibrebound, "reb_build_str").value.decode('ascii')

# Githash
__githash__ = c_char_p.in_dll(clibrebound, "reb_githash_str").value.decode('ascii')

# Check for version
try:
    moduleversion = pkg_resources.require("rebound")[0].version
    libreboundversion = __version__
    if moduleversion != libreboundversion:
        warnings.warn("WARNING: python module and librebound have different version numbers: '%s' vs '%s'.\n" %(moduleversion, libreboundversion), ImportWarning)
except:
    # Might fails on python3 versions, but not important
    pass

# Exceptions    
class SimulationError(Exception):  
    """The simulation exited with a generic error."""
    pass

class Encounter(Exception):
    """The simulation exited because a close encounter has been detected.
    You may want to search for the pair of bodies which have the smallest distance."""
    pass

class Collision(Exception):
    """The simulation exited because a collision has been detected.
    You may want to search for which particles have a lastcollision time equal to the simulation time."""
    pass

class Escape(Exception):
    """The simulation exited because a particle has been se encounter has been detected.
    You may want to search for the particle with the largest distance from the 
    origin and remove it from the simulation."""
    pass

class NoParticles(Exception):
    """The simulation exited because no particles are left in the simulation."""
    pass
    
class ParticleNotFound(Exception):
    """Particle was not found in the simulation."""
    pass

from .simulation import Simulation, Orbit, Variation, reb_simulation_integrator_saba, reb_simulation_integrator_whfast, reb_simulation_integrator_sei, reb_simulation_integrator_mercurius
from .particle import Particle
from .plotting import OrbitPlot
from .tools import hash
from .simulationarchive import SimulationArchive
from .interruptible_pool import InterruptiblePool

__all__ = ["__libpath__", "__version__", "__build__", "__githash__", "SimulationArchive", "Simulation", "Orbit", "OrbitPlot", "Particle", "SimulationError", "Encounter", "Collision", "Escape", "NoParticles", "ParticleNotFound", "InterruptiblePool","Variation", "reb_simulation_integrator_whfast", "reb_simulation_integrator_ias15", "reb_simulation_integrator_saba", "reb_simulation_integrator_sei","reb_simulation_integrator_mercurius", "clibrebound"]
