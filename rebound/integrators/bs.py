from ctypes import c_uint, c_double, c_void_p, CFUNCTYPE, POINTER, c_int, Structure
from .. import clibrebound

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

from ..simulation import Simulation
ODE._fields_ = [
                ("length", c_uint),
                ("y", POINTER(c_double)),
                ("needs_nbody", c_uint),
                ("ref", c_void_p),
                ("_derivatives", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double), c_double)),
                ("_getscale", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double))),
                ("_pre_timestep", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double))),
                ("_post_timestep", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double))),
                ("N_allocated", c_uint),
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

class IntegratorBS(Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_bs.
    It controls the behaviour of the Gragg-Bulirsch-Stoer integrator.
    """
    _fields_ = [
                ("eps_abs", c_double),
                ("eps_rel", c_double),
                ("min_dt", c_double),
                ("max_dt", c_double),
                ("_nbody_ode", POINTER(ODE)),
                ("_sequence", POINTER(c_int)),
                ("_cost_per_step", POINTER(c_int)),
                ("_cost_per_time_unit", POINTER(c_double)),
                ("_optimal_step", POINTER(c_double)),
                ("_coeff", POINTER(c_double)),
                ("dt_proposed", c_double),
                ("_first_or_last_step", c_int),
                ("_previous_rejected", c_int),
                ("_target_iter", c_int),
                ("_user_ode_needs_nbody", c_int),
            ]               
