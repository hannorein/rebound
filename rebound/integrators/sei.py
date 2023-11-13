import ctypes

class IntegratorSEI(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_sei.
    It controls the behaviour of the symplectic SEI integrator for shearing
    sheet calculations. It is described in Rein and Tremaine (2011).
    
    This struct should be accessed via the simulation class only. Here is an 
    example:

    >>> sim = rebound.Simulation()
    >>> sim.ri_sei.OMEGA =  1.58
    
    :ivar float OMEGA:          
        The epicyclic frequency OMEGA. For simulations making use of shearing 
        sheet boundary conditions, REBOUND needs to know the epicyclic frequency. 
        By default, OMEGA is 1. For more details read Rein and Tremaine 2011.
    :ivar float OMEGAZ:          
        The z component of the epicyclic frequency OMEGA. By default, it is assuming
        OMEGAZ is the same as OMEGA.
    """
    _fields_ = [("OMEGA", ctypes.c_double),
                ("OMEGAZ", ctypes.c_double),
                ("_lastdt", ctypes.c_double),
                ("_sindt", ctypes.c_double),
                ("_tandt", ctypes.c_double),
                ("_sindtz", ctypes.c_double),
                ("_tandtz", ctypes.c_double)]

