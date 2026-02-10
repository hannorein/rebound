import ctypes

class IntegratorLeapfrog(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_leapfrog.
    It controls the order of the Leapfrog integrator.

    >>> sim = rebound.Simulation()
    >>> sim.ri_leapfrog.order = 4 # default is 2
    
    :ivar int order
        Sets the order of the integrator. Supported values are 2, 4, 6, and 8.
    """
    _fields_ = [("order", ctypes.c_uint),
                ]

