from ctypes import POINTER, c_int, c_long, c_char_p, byref, pointer
from .simulation import Simulation
from . import clibrebound 
import os
import sys
import math
from collections import Mapping

POINTER_REB_SIM = POINTER(Simulation) 

class SimulationArchive(Mapping):
    """
    SimulationArchive Class.

    This class allows fast access to long running simulations 
    by making use of a SimulationArchive binary file.

    Examples
    --------
    Most simulation parameters can be directly changed with the property syntax:

    >>> sa = rebound.SimulationArchive("archive.bin")
    >>> sim = sa.getSimulation(t=1e6)
    >>> print(sim.particles[1])

    """
    def __getitem__(self, key):
        PY3 = sys.version_info[0] == 3
        if PY3:
            int_types = int,
        else:
            int_types = int, long, 

        if isinstance(key, slice):
            raise AttributeError("Slicing not supported due to optimizations.")
        if not isinstance(key, int_types):
            raise AttributeError("Must access individual simulations with integer index.")
        if key < 0:
            key += len(self)
        if key>= len(self):
            raise IndexError("Index out of range, number of blobs stored in binary: %d."%self.Nblob)
        return self.loadFromBlobAndSynchronize(self, key)

    def loadFromBlobAndSynchronize(self, blob, keep_unsynchronized=1):
        sim = self.simp.contents
        clibrebound.reb_fsr_load_blob.restype = c_int
        fsrlbr = clibrebound.reb_fsr_load_blob(self.simp, self.cfilename, c_long(blob));
        if fsrlbr:
            raise ValueError("Error while loading blob in binary file. Errorcode: %d."%fsrlbr)
        sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
        if self.additional_forces:
            sim.additional_forces = self.additional_forces
        sim.fsr_filename = 0 # Setting this to zero, so no new outputs are generated
        sim.integrator_synchronize()
        return sim

    def __setitem__(self, key, value):
        raise AttributeError("Cannot modify SimulationArchive.")

    def __delitem__(self, key):
        pass

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        return self.Nblob

    def __init__(self,filename,additional_forces=None):
        self.additional_forces = additional_forces
        clibrebound.reb_fsr_restart.restype = POINTER_REB_SIM
        self.cfilename = c_char_p(filename.encode("ascii"))

        # Recreate simulation at t=0
        w = c_int(0)
        clibrebound.reb_create_simulation.restype = POINTER(Simulation)
        simp = clibrebound.reb_create_simulation() 
        clibrebound.reb_create_simulation_from_binary_with_messages(simp, self.cfilename,byref(w))
        if (simp is None) or (w.value & 1):     # Major error
            raise ValueError(BINARY_WARNINGS[0])
        # Note: Other warnings not shown!
        self.simp = simp
        sim = self.simp.contents

        self.filesize = os.path.getsize(filename)
        self.dt = sim.dt
        self.interval = sim.fsr_interval
        self.tmin = 0. # Right now simulations must start at t=0
        self.Nblob = int((self.filesize-sim.fsr_seek_first)/sim.fsr_seek_blob)
        self.tmax = self.tmin + self.interval*(self.Nblob+1)

    def getBlobJustBefore(self, t):
        if t>self.tmax+self.dt or t<self.tmin:
            raise ValueError("Requested time outside of baseline stored in binary fie.")
        bi = max(0,int(math.floor((t-self.dt-self.tmin)/self.interval)))
        bt = self.tmin + self.interval*bi
        return bi, bt

    def getSimulation(self, t, mode='blob', keep_unsynchronized=1):
        """
        Possible values for mode:
         - 'blob' This loads a blob such that sim.t<t.
         - 'close' This integrates the simulation to get to the time t but may overshoot by at most one timestep sim.dt.
         - 'exact' This integrates the simulation to exactly time t. This is not compatible with keep_unsynchronized=1. 
        """
        if mode not in ['blob', 'close', 'exact']:
            raise AttributeError("Unknown mode.")

        bi, bt = self.getBlobJustBefore(t)
        if mode='blob':
            return self.loadFromBlobAndSynchronize(self, bi, keep_unsynchronized=keep_unsynchronized)
        else:
            sim = self.simp.contents
            if sim.t<t and bt-sim.dt<sim.t and sim.ri_whfast.keep_unsynchronized==1:
                # Reuse current simulation
                pass
            else:
                # Load from blob
                clibrebound.reb_fsr_load_blob.restype = c_int
                fsrlbr = clibrebound.reb_fsr_load_blob(self.simp, self.cfilename, c_long(bi));
                if fsrlbr:
                    raise ValueError("Error while loading blob in binary file. Errorcode: %d."%fsrlbr)

            if mode=='exact':
                keep_unsynchronized==0
            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            if self.additional_forces:
                sim.additional_forces = self.additional_forces
            sim.fsr_filename = 0 # Setting this to zero, so no new outputs are generated
            exact_finish_time = 1 if mode=='exact' else 0
            sim.integrate(t,exact_finish_time=exact_finish_time)
                
            return sim

    def getSimulations(self, times, mode='blob', keep_unsynchronized=1):
        times.sort()
        for t in times:
            yield self.getSimulation(t, mode=mode, keep_unsynchronized=keep_unsynchronized)

