import sys
from ctypes import c_uint32, c_uint, c_uint64, addressof, POINTER, pointer, byref
from .particle import Particle
from . import clibrebound, ParticleNotFound

try:
    # Required for Python>=3.9
    from collections.abc import MutableMapping
except:
    from collections import MutableMapping


class Particles(MutableMapping):
    """
    This class allows the user to access particles like a dictionary using the particle's index or name.
    Allows for negative indices and slicing.
    """
    def __init__(self, sim):
        self.sim = sim

    @property
    def _ps(self):
        ParticleList = Particle*self.sim.N
        pl = ParticleList.from_address(addressof(self.sim._particles.contents))
        pl._sim = self.sim # keep reference to sim until ParticleList is deallocated to avoid memory issues
        return pl

    def __getitem__(self, key):
        string_types = str,
        int_types = int,
        try:
            import numpy as np
            int_types += np.int64,
        except:
            pass
       
        if isinstance(key, slice):
            return [self[i] for i in range(*key.indices(len(self)))]

        if isinstance(key, int_types):
            if key < 0: # accept negative indices
                key += self.sim.N
            if key < 0 or key >= self.sim.N:
                raise AttributeError("Index {0} used to access particles out of range.".format(key))
            return self._ps[key]

        if isinstance(key, string_types):
            clibrebound.reb_simulation_get_particle_by_name.restype = POINTER(Particle)
            s = key.encode("utf-8")
            ptr = clibrebound.reb_simulation_get_particle_by_name(byref(self.sim), s) 
            
            if ptr:
                p = Particle
                return p.from_address(addressof(ptr.contents))
            else:
                raise ParticleNotFound("Particle with name \"{0}\" was not found in the simulation.".format(key)) 
        raise AttributeError("Unable to get particles with key {0}.".format(key))

    def __setitem__(self, key, value):
        if isinstance(value, Particle):
            value._sim = pointer(self.sim)
            p = self[key]
            if p.index == -1:
                raise AttributeError("Can't set particle (particle not found in simulation).")
            else:
                self._ps[p.index] = value

    def __delitem__(self, key):
        pass

    def __iter__(self):
        if self.sim.N>0:
            for p in self._ps:
                yield p

    def __len__(self):
        return self.sim.N
