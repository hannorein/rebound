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
        if hasattr(self, "_needsfree"): # Root node needs to free its children
            clibrebound.reb_orbit_hierarchy_free(byref(self))

    def print(self, sim):
        clibrebound.reb_orbit_hierarchy_print(byref(self), byref(sim), c_int(0))

        
    def is_jacobi(self):
        clibrebound.reb_orbit_hierarchy_is_jacobi.restype = c_int
        return bool(clibrebound.reb_orbit_hierarchy_is_jacobi(byref(self)))
    
    def is_jacobi_ordered(self, sim):
        clibrebound.reb_orbit_hierarchy_is_jacobi_ordered.restype = c_int
        return bool(clibrebound.reb_orbit_hierarchy_is_jacobi_ordered(byref(self), byref(sim)))

    def orbit_generator(self):
        if self.primary and self.secondary:
            p1 = self.secondary.contents.com.contents
            p1m = p1.m
            p2 = self.primary.contents.com.contents
            p2m = p2.m
            com = self.com.contents
            com.vx = 0.0
            com.vy = 0.0
            com.vz = 0.0
            p = p1-p2+com
            p.m = (p1m*p2m)/(p1m+p2m)
            fac = p2m/(p1m+p2m)
            com.m = p1m+p2m-p.m
            yield p, com, fac
            p = p2-p1+com
            p.m = (p1m*p2m)/(p1m+p2m)
            fac = p1m/(p1m+p2m)
            com.m = p1m+p2m-p.m
            yield p, com, fac

        if self.primary:
            yield from self.primary.contents.orbit_generator()
        if self.secondary:
            yield from self.secondary.contents.orbit_generator()
    
    def particle_generator(self):
        if self.primary and self.secondary:
            yield from self.primary.contents.particle_generator()
            yield from self.secondary.contents.particle_generator()
        else:
            yield self.com.contents


    def Norbits(self):
        n = 0
        if self.primary:
            n = n + 1 + self.primary.contents.Norbits()
        if self.secondary:
            n = n + 1 + self.secondary.contents.Norbits()
        return n

OrbitHierarchy._fields_ = [
        ("primary", POINTER(OrbitHierarchy)),
        ("secondary", POINTER(OrbitHierarchy)),
        ("com", POINTER(Particle)),
        ]

