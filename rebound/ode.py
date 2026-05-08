from . import clibrebound
from ctypes import Structure, c_double, POINTER, c_uint, CFUNCTYPE, c_size_t, c_void_p 
class ODE(Structure):
    @property
    def derivatives(self):
        raise AttributeError("You can only set C function pointers from python.")
    @derivatives.setter
    def derivatives(self, func):
        self._dfp = ODEDER(func)
        func.argtypes = self._dfp.argtypes # I do not understand why this is needed
        self._derivatives = self._dfp
    def update_particles(self):
        clibrebound.reb_integrator_bs_update_particles(self.r, None) 

from .simulation import Simulation
ODE._fields_ = [
                ("length", c_uint),
                ("y", POINTER(c_double)),
                ("needs_nbody", c_uint),
                ("ref", c_void_p),
                ("_derivatives", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double), c_double)),
                ("_getscale", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double))),
                ("_pre_timestep", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double))),
                ("_post_timestep", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double))),
                ("N_allocated", c_size_t),
                ("_scale", POINTER(c_double)),
                ("_C", POINTER(c_double)),
                ("_D", POINTER(POINTER(c_double))),
                ("_y1", POINTER(c_double)),
                ("_y0Dot", POINTER(c_double)),
                ("_yDot", POINTER(c_double)),
                ("_yTmp", POINTER(c_double)),
                ("r", POINTER(Simulation)),
            ]               

ODEDER = CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double), c_double)
