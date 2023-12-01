
import ctypes

from .. import clibrebound
from ..simulation import Simulation
from ..particle import Particle
from ..vectors import Vec3dBasic
    
class IntegratorTRACE(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_trace.
    It controls the behaviour of the TRACE integrator.  See Lu et al. (2023)
    for more details.

    FIX THIS
    :ivar float hillfac:
        Switching radius in units of the hill radius.

    Example usage:

    >>> sim = rebound.Simulation()
    >>> sim.integrator = "trace"
    >>> sim.ri_trace.hillfac = 4.

    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, is_synchronized={4}, hillfac={5}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.is_synchronized, self.hillfac)

    _fields_ = [("_S", ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(Simulation), ctypes.c_uint, ctypes.c_uint)),
                ("_S_peri", ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(Simulation), ctypes.c_uint, ctypes.c_uint)),
                ("hillfac", ctypes.c_double),
                ("peri", ctypes.c_double),
                ("vfac_p", ctypes.c_double),
                ("recalculate_coordinates_this_timestep", ctypes.c_uint),
                ("safe_mode", ctypes.c_uint),
                ("is_synchronized", ctypes.c_uint),
                ("mode", ctypes.c_uint),
                ("_encounterN", ctypes.c_uint),
                ("_encounterNactive", ctypes.c_uint),
                ("_allocated_N", ctypes.c_uint),
                ("_allocated_N_additionalforces", ctypes.c_uint),
                ("_particles_backup", ctypes.POINTER(Particle)),
                ("_particles_backup_try", ctypes.POINTER(Particle)),
                ("_encounter_map", ctypes.POINTER(ctypes.c_int)),
                ("_encounter_map_internal", ctypes.POINTER(ctypes.c_int)),
                ("_com_pos", Vec3dBasic),
                ("_com_vel", Vec3dBasic),
                ("_current_Ks", ctypes.POINTER(ctypes.c_int)),
                ("current_L", ctypes.c_uint),
                ("print", ctypes.c_uint),
                ]
    @property
    def S(self):
        raise AttributeError("You can only set C function pointers from python.")
    @S.setter
    def S(self, func):
        if func == "default":
            self._S = cast(clibrebound.reb_integrator_trace_switch_default,TRACEKF)
        else:
            self._Sfp = TRACEKF(func) # what is this
            self._S = self._Sfp

    @property
    def S_peri(self):
        raise AttributeError("You can only set C function pointers from python.")
    @S_peri.setter
    def S_peri(self, func):
        if func == "default":
            self._S_peri = cast(clibrebound.reb_integrator_trace_switch_fdot_peri,TRACELF)
        elif func == "distance":
            self._S_peri = cast(clibrebound.reb_integrator_trace_peri_switch_default,TRACELF)
        else:
            self._S_perifp = TRACELF(func) # what is this
            self._S_peri = self._S_perifp

TRACEKF = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(Simulation), ctypes.c_double, ctypes.c_double)
TRACELF = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(Simulation), ctypes.c_double)

