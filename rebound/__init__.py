# -*- coding: utf-8 -*-
"""An N-body integrator package for python."""

import sys
import os
import warnings
import platform
from ctypes import cdll, c_char_p

# Find suffix
if platform.system()=="Windows":
    # Using distutils.sysconfig instead of sysconfig because 
    # of a bug in Python < 3.8 on windows
    import distutils.sysconfig as sysconfig
else:
    import sysconfig
suffix = sysconfig.get_config_var('EXT_SUFFIX')

if suffix is None:
    suffix = ".so"

try: # Only needed for pyodide
    import pyodide_js
    from site import getsitepackages
    pyodide_js._module.loadDynamicLibrary(getsitepackages()[0]+"/librebound"+suffix)
    del getsitepackages
    del pyodide_js
except:
    pass


# Make changes for python 2 and 3 compatibility
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x


# Import shared library
pymodulepath = os.path.dirname(os.path.abspath(__file__))
pymodulepath = os.path.abspath(os.path.join(pymodulepath, os.pardir))
__libpath__ = os.path.join(pymodulepath, "librebound"+suffix)
clibrebound = cdll.LoadLibrary(__libpath__)

# Version
__version__ = c_char_p.in_dll(clibrebound, "reb_version_str").value.decode('ascii')

# Build
__build__ = c_char_p.in_dll(clibrebound, "reb_build_str").value.decode('ascii')

# Githash
__githash__ = c_char_p.in_dll(clibrebound, "reb_githash_str").value.decode('ascii')

# Check for version
try:
    import pkg_resources
    moduleversion = pkg_resources.require("rebound")[0].version
    libreboundversion = __version__
    if moduleversion != libreboundversion:
        warnings.warn("WARNING: python module and librebound have different version numbers: '%s' vs '%s'.\n" %(moduleversion, libreboundversion), ImportWarning)
except:
    # Might fail in some python3 setups, but not important
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

from .tools import hash, mod2pi, M_to_f, E_to_f, M_to_E, spherical_to_xyz, xyz_to_spherical
from .simulation import Simulation, Orbit, Variation, reb_simulation_integrator_saba, reb_simulation_integrator_whfast, reb_simulation_integrator_sei, reb_simulation_integrator_mercurius, reb_simulation_integrator_ias15, ODE, Rotation, Vec3d, _Vec3d, binary_field_descriptor_list
from .particle import Particle
from .plotting import OrbitPlot, OrbitPlotSet
from .simulationarchive import SimulationArchive

class InterruptiblePool:
    def __init__(self, processes=None, initializer=None, initargs=(), **kwargs):
        raise RuntimeError("InterruptiblePool is deprecated. Use the multiprocess module instead.")

__all__ = ["__libpath__", "__version__", "__build__", "__githash__", "SimulationArchive", "Simulation", "Orbit", "OrbitPlot", "OrbitPlotSet", "Particle", "SimulationError", "Encounter", "Collision", "Escape", "NoParticles", "ParticleNotFound", "InterruptiblePool","Variation", "reb_simulation_integrator_whfast", "reb_simulation_integrator_ias15", "reb_simulation_integrator_saba", "reb_simulation_integrator_sei","reb_simulation_integrator_mercurius", "clibrebound", "mod2pi", "M_to_f", "E_to_f", "M_to_E", "ODE", "Rotation", "Vec3d", "spherical_to_xyz", "xyz_to_spherical","binary_field_descriptor_list"]
