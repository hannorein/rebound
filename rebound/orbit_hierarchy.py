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

    def print(self, sim, level=0):
        space = "  "*level
        if self.primary and self.secondary:
            print(space+"binary m=%e"%self.com.contents.m)
            self.primary.contents.print(sim,level=level+1)
            self.secondary.contents.print(sim,level=level+1)
        else:
            print(space+"particle m=%e"%self.com.contents.m)

        
    def is_jacobi(self):
        clibrebound.reb_orbit_hierarchy_is_jacobi.restype = c_int
        return bool(clibrebound.reb_orbit_hierarchy_is_jacobi(byref(self)))
    
    def is_jacobi_ordered(self, sim):
        clibrebound.reb_orbit_hierarchy_is_jacobi_ordered.restype = c_int
        return bool(clibrebound.reb_orbit_hierarchy_is_jacobi_ordered(byref(self), byref(sim)))

    def orbit_generator(self, style="complex",star=None):
        if style=="complex":
            if self.primary and self.secondary:
                p1m = self.primary.contents.com.contents.m
                p2m = self.secondary.contents.com.contents.m
                com = self.com.contents
                G = com._sim.contents.G
                p1 = self.primary.contents.com.contents 
                p2 = self.secondary.contents.com.contents
                o1 = p1.orbit(primary=p2, G=G)
                o1.a *= p2m/(p1m+p2m) 
                o2 = p2.orbit(primary=p1, G=G)
                o2.a *= p1m/(p1m+p2m) 
                yield o1, com, p1
                yield o2, com, p2

            if self.primary:
                yield from self.primary.contents.orbit_generator(style=style)
            if self.secondary:
                yield from self.secondary.contents.orbit_generator(style=style)
        if style=="jacobi":
            if not self.is_jacobi():
                raise RuntimeError("Cannot generate Jacobi orbits as is_jacobi()==False.")
            if self.primary and self.secondary:
                p1m = self.primary.contents.com.contents.m
                p2m = self.secondary.contents.com.contents.m
                G = self.com.contents._sim.contents.G
                p1 = self.primary.contents.com.contents 
                p2 = self.secondary.contents.com.contents
                o2 = p2.orbit(primary=p1, G=G)
                yield from self.primary.contents.orbit_generator(style=style)
                yield o2, p1, p2
        if style=="heliocentric":
            if star==None:
                star = self.com.contents._sim.contents.particles[0]
            if self.primary and self.secondary:
                yield from self.primary.contents.orbit_generator(style=style,star=star)
                yield from self.secondary.contents.orbit_generator(style=style,star=star)
            else:
                G = self.com.contents._sim.contents.G
                if star != self.com.contents:
                    o = self.com.contents.orbit(primary=star, G=G)
                    yield o, star, self.com.contents

    
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

