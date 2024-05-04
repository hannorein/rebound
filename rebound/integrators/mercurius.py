import ctypes

from .. import clibrebound
from ..simulation import Simulation
from ..particle import Particle
from ..vectors import Vec3dBasic

class IntegratorMercurius(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_mercurius.
    It controls the behaviour of the MERCURIUS integrator.  See Rein et al. (2019) 
    for more details.
    
    :ivar float r_crit_hill:      
        Switching radius in units of the hill radius.

    Example usage:
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "mercurius"
    >>> sim.ri_mercurius.r_crit_hill = 3.

    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, is_synchronized={4}, r_crit_hill={5}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.is_synchronized, self.r_crit_hill)

    _fields_ = [("_L", ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(Simulation), ctypes.c_double, ctypes.c_double)),
                ("r_crit_hill", ctypes.c_double),
                ("recalculate_coordinates_this_timestep", ctypes.c_uint),
                ("recalculate_r_crit_this_timestep", ctypes.c_uint),
                ("safe_mode", ctypes.c_uint),
                ("is_synchronized", ctypes.c_uint),
                ("mode", ctypes.c_uint),
                ("_encounter_N", ctypes.c_uint),
                ("_encounter_N_active", ctypes.c_uint),
                ("_tponly_encounter", ctypes.c_uint),
                ("_N_allocated", ctypes.c_uint),
                ("_N_allocated_additional_forces", ctypes.c_uint),
                ("_N_allocated_dcrit", ctypes.c_uint),
                ("_dcrit", ctypes.POINTER(ctypes.c_double)),
                ("_particles_backup", ctypes.POINTER(Particle)),
                ("_particles_backup_additional_forces", ctypes.POINTER(Particle)),
                ("_encounter_map", ctypes.POINTER(ctypes.c_int)),
                ("_com_pos", Vec3dBasic),
                ("_com_vel", Vec3dBasic),
                ]
    @property
    def L(self):
        raise AttributeError("You can only set C function pointers from python.")
    @L.setter
    def L(self, func):
        if func == "mercury":
            self._L = ctypes.cast(clibrebound.reb_integrator_mercurius_L_mercury,MERCURIUSLF)
        elif func == "C4":
            self._L = ctypes.cast(clibrebound.reb_integrator_mercurius_L_C4,MERCURIUSLF)
        elif func == "C5":
            self._L = ctypes.cast(clibrebound.reb_integrator_mercurius_L_C5,MERCURIUSLF)
        elif func == "infinity":
            self._L = ctypes.cast(clibrebound.reb_integrator_mercurius_L_infinity,MERCURIUSLF)
        else:
            self._Lfp = MERCURIUSLF(func)
            self._L = self._Lfp

MERCURIUSLF = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.POINTER(Simulation), ctypes.c_double, ctypes.c_double)
