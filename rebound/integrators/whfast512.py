import ctypes
from ..particle import Particle

WHFAST512_COORDINATES = {"jacobi": 0, "democraticheliocentric": 1}

class IntegratorWHFast512(ctypes.Structure):
    _fields_ = [("gr_potential", ctypes.c_uint),
                ("N_systems", ctypes.c_uint),
                ("keep_unsynchronized", ctypes.c_uint),
                ("concatenate_steps", ctypes.c_uint),
                ("corrector", ctypes.c_uint),
                ("_coordinates", ctypes.c_int),
                ("_is_synchronized", ctypes.c_uint),
                ("_N_allocated_particles_keep_unsynchronized", ctypes.c_uint),
                ("_particles_keep_unsynchronized", ctypes.POINTER(Particle)),
                ("_recalculate_constants", ctypes.c_uint),
                ("_time_of_last_synchronize", ctypes.c_double),
                ("_N_allocated", ctypes.c_uint),
                ("_p512", ctypes.POINTER(Particle)),  # Note: actually a pointer to reb_particle_avx512
                ("mat8_0", ctypes.POINTER(ctypes.c_double)),
                ("mat8_1", ctypes.POINTER(ctypes.c_double)),
                ("mat8_2", ctypes.POINTER(ctypes.c_double)),
                ]
    
    @property
    def coordinates(self):
        """
        Get or set the internal coordinate system.

        Available coordinate systems are:

        - ``'jacobi'`` (default)
        - ``'democraticheliocentric'``
        """
        i = self._coordinates
        for name, _i in WHFAST512_COORDINATES.items():
            if i==_i:
                return name
        return i
    @coordinates.setter
    def coordinates(self, value):
        if isinstance(value, int):
            self._coordinates = ctypes.c_uint(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in WHFAST512_COORDINATES: 
                self._coordinates = WHFAST512_COORDINATES[value]
            else:
                raise ValueError("Warning. Coordinate system not found.")
