import sys
from ctypes import c_uint32, c_uint, c_uint64, addressof, POINTER, pointer, byref
from .hash import hash as rebhash
from .particle import Particle
from . import clibrebound, ParticleNotFound

try:
    # Required for Python>=3.9
    from collections.abc import MutableMapping
except:
    from collections import MutableMapping


class Particles(MutableMapping):
    """
    This class allows the user to access particles like a dictionary using the particle's 1) index 2) hash 3) string (which will be converted to hash).
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
        hash_types = c_uint32, c_uint, c_uint64
        PY3 = sys.version_info[0] == 3
        if PY3:
            string_types = str,
            int_types = int,
        else:
            string_types = basestring,
            int_types = int, long,
       
        if isinstance(key, slice):
            return [self[i] for i in range(*key.indices(len(self)))]

        if isinstance(key, int_types):
            if key < 0: # accept negative indices
                key += self.sim.N
            if key < 0 or key >= self.sim.N:
                raise AttributeError("Index {0} used to access particles out of range.".format(key))
            return self._ps[key]

        else:
            clibrebound.reb_simulation_particle_by_hash.restype = POINTER(Particle)
            if isinstance(key, string_types):
                key = rebhash(key)
            elif not isinstance(key, hash_types):
                raise AttributeError("Expecting string, integer or ctypes.c_uint32 as argument to sim.particles.  See UniquelyIdentifyingParticlesWithHashes.ipynb ipython_example.")

            ptr = clibrebound.reb_simulation_particle_by_hash(byref(self.sim), key) 
            
            if ptr:
                p = Particle
                return p.from_address(addressof(ptr.contents))
            else:
                raise ParticleNotFound("Particle was not found in the simulation.") 

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
