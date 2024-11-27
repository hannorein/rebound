import ctypes
from ..particle import Particle

WHFAST_KERNELS = {"default": 0, "modifiedkick": 1, "composition": 2, "lazy": 3}
WHFAST_COORDINATES = {"jacobi": 0, "democraticheliocentric": 1, "whds": 2, "barycentric": 3}

class IntegratorWHFast(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_whfast.
    It controls the behaviour of the symplectic WHFast integrator described 
    in Rein and Tamayo (2015) and in Rein, Tamayo, and Brown (2019).
    
    This struct should be accessed via the simulation class only. Here is an 
    example:

    >>> sim = rebound.Simulation()
    >>> sim.integrator = "whfast"
    >>> sim.ri_whfast.corrector =  11
    >>> sim.ri_whfast.kernel = "lazy"

    
    :ivar int corrector:      
        The order of the symplectic corrector in the WHFast integrator.
        By default, the symplectic correctors are turned off (=0). For high
        accuracy simulation set this value to 11 or 17. For more details read 
        Rein and Tamayo (2015).
    :ivar int corrector2:      
        Second correctors (C2 of Wisdom et al 1996).
        By default, the second symplectic correctors are turned off (=0). 
        Set to 1 to turn them on.
    :ivar int/string kernel:      
        Kernel option. Set to 0 for the default WH kernel (standard kick step).
        Other options are "modifiedkick" (1), "composition" (2), "lazy" (3).
    :ivar int recalculate_coordinates_this_timestep:
        Sets a flag that tells WHFast that the particles have changed.
        Setting this flag to 1 (default 0) triggers the WHFast integrator
        to recalculate Jacobi/heliocenctric coordinates. This is needed 
        if the user changes the particle position, velocity or mass 
        in-between timesteps.  After every timestep the flag is set back 
        to 0, so if you continuously update the particles manually, 
        you need to set this flag to 1 after every timestep.
    :ivar string coordinates:
        Sets the internal coordinate system that WHFast is using. By default
        it uses ``'jacobi'`` (=0) coordinates. Other options are 
        ``'democraticheliocentric'`` (=1) and ``'whds'`` (=2). See Hernandez 
        and Dehnen (2017) for more information.
    :ivar int safe_mode:
        If safe_mode is 1 (default) particles can be modified between
        timesteps and particle velocities and positions are always synchronised.
        If you set safe_mode to 0, the speed and accuracy of WHFast improve.
        However, make sure you are aware of the consequences. Read the iPython tutorial
        on advanced WHFast usage to learn more.
    """
    _fields_ = [("corrector", ctypes.c_uint),
                ("corrector2", ctypes.c_uint),
                ("_kernel", ctypes.c_uint),
                ("_coordinates", ctypes.c_uint),
                ("recalculate_coordinates_this_timestep", ctypes.c_uint),
                ("safe_mode", ctypes.c_uint),
                ("keep_unsynchronized", ctypes.c_uint),
                ("_p_jh", ctypes.POINTER(Particle)),
                ("_p_temp", ctypes.POINTER(Particle)),
                ("is_synchronized", ctypes.c_uint),
                ("_N_allocated", ctypes.c_uint),
                ("_N_allocated_tmp", ctypes.c_uint),
                ("_timestep_warning", ctypes.c_uint),
                ("_recalculate_coordinates_but_not_synchronized_warning", ctypes.c_uint)]

    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, keep_unsynchonized={4}, is_synchronized={5}, corrector={6}, corrector2={7}, kernel={8}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.keep_unsynchronized, self.is_synchronized, self.corrector, self.corrector2, self.kernel)

    @property
    def coordinates(self):
        """
        Get or set the internal coordinate system.

        Available coordinate systems are:

        - ``'jacobi'`` (default)
        - ``'democraticheliocentric'``
        - ``'whds'``
        - ``'barycentric'``
        """
        i = self._coordinates
        for name, _i in WHFAST_COORDINATES.items():
            if i==_i:
                return name
        return i
    @coordinates.setter
    def coordinates(self, value):
        if isinstance(value, int):
            self._coordinates = ctypes.c_uint(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in WHFAST_COORDINATES: 
                self._coordinates = WHFAST_COORDINATES[value]
            else:
                raise ValueError("Warning. Coordinate system not found.")
    @property
    def kernel(self):
        """
        Get or set the WHFast Kernel.

        Available kernels are:

        - ``'default'`` (standard WH kernel, kick)
        - ``'modifiedkick'`` (modified kick for newtonian gravity)
        - ``'composition'`` (Wisdom's composition method)
        - ``'lazy'`` (Lazy implementer's method)
        """
        i = self._kernel
        for name, _i in WHFAST_KERNELS.items():
            if i==_i:
                return name
        return i
    @kernel.setter
    def kernel(self, value):
        if isinstance(value, int):
            self._kernel = ctypes.c_uint(value)
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "")
            if value in WHFAST_KERNELS: 
                self._kernel = WHFAST_KERNELS[value]
            else:
                raise ValueError("Warning. Kernel not found.")


