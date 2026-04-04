import ctypes

class IntegratorCustom(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator.
    It controls the behaviour of a custom, user-provided integrator.
    """
    _fields_ = [("id", ctypes.c_uint32),
                ("step", ctypes.c_void_p),
                ("synchronize", ctypes.c_void_p),
                ("reset", ctypes.c_void_p),
                ("data", ctypes.c_void_p),
                ("data_size", ctypes.c_size_t),
                ]

