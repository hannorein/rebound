import ctypes

from ..particle import Particle

class IntegratorWHFast512(ctypes.Structure):
    _fields_ = [("gr_potential", ctypes.c_uint),
                ("N_systems", ctypes.c_uint),
                ("keep_unsynchronized", ctypes.c_uint),
                ("is_synchronized", ctypes.c_uint),
                ("_N_allocated", ctypes.c_uint),
                ("recalculate_constants", ctypes.c_uint),
                ("_dt_cached", ctypes.c_double),
                ("coordinates", ctypes.c_int),
                ("_p512", ctypes.POINTER(Particle)),
                ("_p_jh0", Particle*4),
                ("mat8_0", ctypes.POINTER(ctypes.c_double)),
                ("mat8_1", ctypes.POINTER(ctypes.c_double)),
                ("mat8_2", ctypes.POINTER(ctypes.c_double)),
                ]
