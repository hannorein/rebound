from ctypes import POINTER, Structure, c_int, byref
from . import clibrebound 
from .particle import Particle

class OrbitHierarchy(Structure):
    """
    A class representing one orbit in a tree. 
    """
    def __repr__(self):
        return "<rebound.OrbitHierarchy instance, primary={0} secondary={1} com={2}>".format(str(self.primary),str(self.secondary), str(self.com))
    
    def __del__(self):
        if self._b_needsfree_ == 0: # Allocated from within C
            clibrebound.reb_orbit_hierarchy_free(byref(self))

    def print(self, sim):
        clibrebound.reb_orbit_hierarchy_print(byref(self), byref(sim), c_int(0))

        
    def is_jacobi(self):
        clibrebound.reb_orbit_hierarchy_is_jacobi.restype = c_int
        return bool(clibrebound.reb_orbit_hierarchy_is_jacobi(byref(self)))
    
    def is_jacobi_ordered(self, sim):
        clibrebound.reb_orbit_hierarchy_is_jacobi_ordered.restype = c_int
        return bool(clibrebound.reb_orbit_hierarchy_is_jacobi_ordered(byref(self), byref(sim)))

OrbitHierarchy._fields_ = [
        ("primary", POINTER(OrbitHierarchy)),
        ("secondary", POINTER(OrbitHierarchy)),
        ("com", POINTER(Particle)),
        ]

