import ctypes

from ..particle import Particle

class IntegratorWHFast512(ctypes.Structure):
    _fields_ = [("gr_potential", ctypes.c_uint),
                ("N_systems", ctypes.c_uint),
                ("keep_unsynchronized", ctypes.c_uint),
                ("is_synchronized", ctypes.c_uint),
                ("_N_allocated", ctypes.c_uint),
                ("recalculate_constants", ctypes.c_uint),
                ("_p_jh", ctypes.POINTER(Particle)),
                ("_p_jh0", Particle*4)]
