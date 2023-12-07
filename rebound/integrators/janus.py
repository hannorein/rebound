import ctypes

class ParticleInt(ctypes.Structure):
    """ Used for Janus integrator only """
    _fields_ = [
                ("x", ctypes.c_int64),
                ("y", ctypes.c_int64),
                ("z", ctypes.c_int64),
                ("vx", ctypes.c_int64),
                ("vy", ctypes.c_int64),
                ("vz", ctypes.c_int64),
                ]

class IntegratorJanus(ctypes.Structure):
    _fields_ = [
                ("scale_pos",ctypes.c_double),
                ("scale_vel",ctypes.c_double),
                ("order", ctypes.c_uint),
                ("recalculate_integer_coordinates_this_timestep", ctypes.c_uint),
                ("p_int", ctypes.POINTER(ParticleInt)),
                ("_N_allocated",ctypes.c_uint),
                ]
