import ctypes
from ctypes import cast

from .. import clibrebound
from ..simulation import Simulation
from ..particle import Particle
from ..vectors import Vec3dBasic

BRACE_PERI_MODES = {"FULL_BS": 1, "PARTIAL_BS": 0, "FULL_IAS15": 2}

class IntegratorBRACE(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_brace.
    It controls the behaviour of the BRACE integrator.  See Lu, Hernandez and Rein (2024)
    for more details.

    :ivar float r_crit_hill:
        Switching radius in units of the modified hill radius.
    :ivar float peri_mode:
        Pericenter switching mode. Determines how close encounters with the central star are integrated.
        The default is 1, which uses the FULL BS prescription. This converts our system from Democratic Heliocentric Coordinates to back to Intertial Cartesian coordinates. The whole system is integrated with Bulirsch-Stoer.
        If set to 0, we use the PARTIAL BS prescription. This moves all terms from the Jump step to the Kepler step as described in Hernandez and Dehnen (2023). The Kepler step is integrated with Bulirsch-Stoer.
        If set to 2, we use the FULL IAS15 prescription. This acts similarly to the PURE BS prescription, but uses IAS15 instead of BS to integrate the whole system.
    :ivar float peri_crit_eta:
        Pericenter switching condition: ratio of timestep to condition from Pham, Rein and Spiegel 2024.

    Example usage:

    >>> sim = rebound.Simulation()
    >>> sim.integrator = "brace"
    >>> sim.ri_brace.r_crit_hill = 4
    >>> sim.ri_brace.peri_crit_eta = 1

    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, r_crit_hill={3}, peri_mode=={4}, peri_crit_eta=={5}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.r_crit_hill, self.peri_mode, self.peri_crit_eta)

    _fields_ = [("_S", ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(Simulation), ctypes.c_uint, ctypes.c_uint)),
                ("_S_peri", ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(Simulation), ctypes.c_uint)),
                ("r_crit_hill", ctypes.c_double),
                ("peri_crit_eta", ctypes.c_double),
                ("old_t", ctypes.c_double),
                ("nextias_dt", ctypes.c_double),
                ("_mode", ctypes.c_uint),
                ("_encounter_integrator", ctypes.c_uint),
                ("_encounter_N", ctypes.c_uint),
                ("_encounter_N_active", ctypes.c_uint),
                ("_N_allocated", ctypes.c_uint),
                ("_N_allocated_additionalforces", ctypes.c_uint),
                ("_particles_backup", ctypes.POINTER(Particle)),
                ("_particles_backup_additional_forces", ctypes.POINTER(Particle)),
                ("_encounter_map", ctypes.POINTER(ctypes.c_int)),
                ("_com_pos", Vec3dBasic),
                ("_com_vel", Vec3dBasic),
                ("_current_Ks", ctypes.POINTER(ctypes.c_int)),
                ("_force_accept", ctypes.c_uint),
                ]
    @property
    def S(self):
        raise AttributeError("You can only set C function pointers from python.")
    @S.setter
    def S(self, func):
        if func == "default":
            self._S = cast(clibrebound.reb_integrator_brace_switch_default,BRACEKF)
        else:
            self._Sfp = BRACEKF(func)
            self._S = self._Sfp

    @property
    def S_peri(self):
        raise AttributeError("You can only set C function pointers from python.")
    @S_peri.setter
    def S_peri(self, func):
        if func == "default":
            self._S_peri = cast(clibrebound.reb_integrator_brace_switch_peri_default,BRACECF)
        elif func == "none":
            self._S_peri = cast(clibrebound.reb_integrator_brace_switch_peri_none,BRACECF)
        else:
            self._S_perifp = BRACECF(func)
            self._S_peri = self._S_perifp
    
    @property
    def peri_mode(self):
        """
        Get or set the pericenter switching mode.

        Available kernels are:

        - ``'FULL_BS'`` (Integrate entire system with BS)
        - ``'PARTIAL_BS'`` (Integrate only the Kepler step with BS)
        - ``'FULL_IAS15'`` (Integrate entire system with IAS15)
        """
        i = self._peri_mode
        for name, _i in BRACE_PERI_MODES.items():
            if i==_i:
                return name
        return i
    @peri_mode.setter
    def peri_mode(self, value):
        if isinstance(value, int):
            self._peri_mode = ctypes.c_uint(value)
        elif isinstance(value, basestring):
            if value in BRACE_PERI_MODES: 
                self._peri_mode = BRACE_PERI_MODES[value]
            else:
                raise ValueError("Warning. Pericenter switching mode not found.")

BRACEKF = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(Simulation), ctypes.c_uint, ctypes.c_uint)
BRACECF = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(Simulation), ctypes.c_uint)
