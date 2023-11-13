import ctypes


class reb_dp7(ctypes.Structure):
    _fields_ = [("p0", ctypes.POINTER(ctypes.c_double)),
                ("p1", ctypes.POINTER(ctypes.c_double)),
                ("p2", ctypes.POINTER(ctypes.c_double)),
                ("p3", ctypes.POINTER(ctypes.c_double)),
                ("p4", ctypes.POINTER(ctypes.c_double)),
                ("p5", ctypes.POINTER(ctypes.c_double)),
                ("p6", ctypes.POINTER(ctypes.c_double))]

class IntegratorIAS15(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_ias15.
    It controls the behaviour of the SEI integrator. See Rein & Spiegel (2015)
    for more information.

    Example usage:

    >>> sim = rebound.Simulation()
    >>> sim.integrator = "ias15"
    >>> sim.ri_ias15.epsilon = 0.

    :ivar float epsilon:          
        Controls the precision of the integrator. Set to 0 for fixed timesteps.
    
    :ivar float min_dt:          
        IAS15 is an adaptive method. This sets the minimum timestep.
    
    :ivar float adaptive_mode:
        Determines how the adaptive timestep is chosen. 
        This replaces the previous epsilon_global variable.
        TODO: list options.
    
    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, epsilon={3}, min_dt={4}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.epsilon, self.min_dt)
    
    _fields_ = [("epsilon", ctypes.c_double),
                ("min_dt", ctypes.c_double),
                ("adaptive_mode", ctypes.c_uint),
                ("_iterations_max_exceeded", ctypes.c_uint64),
                ("_N_allocated", ctypes.c_uint),
                ("_at", ctypes.POINTER(ctypes.c_double)),
                ("_x0", ctypes.POINTER(ctypes.c_double)),
                ("_v0", ctypes.POINTER(ctypes.c_double)),
                ("_a0", ctypes.POINTER(ctypes.c_double)),
                ("_csx", ctypes.POINTER(ctypes.c_double)),
                ("_csv", ctypes.POINTER(ctypes.c_double)),
                ("_csa0", ctypes.POINTER(ctypes.c_double)),
                ("_g", reb_dp7),
                ("_b", reb_dp7),
                ("_csb", reb_dp7),
                ("_e", reb_dp7),
                ("_br", reb_dp7),
                ("_er", reb_dp7),
                ("_map", ctypes.POINTER(ctypes.c_int)),
                ("_map_allocated_n", ctypes.c_uint),
                ]
