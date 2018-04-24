from ctypes import Structure, c_double, POINTER, c_float, c_int, c_uint, c_uint32, c_int64, c_long, c_ulong, c_ulonglong, c_void_p, c_char_p, CFUNCTYPE, byref, create_string_buffer, addressof, pointer, cast
from .simulation import Simulation, BINARY_WARNINGS
from . import clibrebound 
import os
import sys
import math
import warnings

POINTER_REB_SIM = POINTER(Simulation) 

class SimulationArchive(Structure):
    """
    SimulationArchive Class.

    The SimulationArchive is a binary file that includes all
    settings, constants as well as initial conditions to reproduce
    a simulation exactly (down to the last bit). It is also possible 
    to add snapshots of the current state as the simulation runs.
    This allows for fast access to long running simulations 
    by making use of the SimulationArchive binary file. For a full
    discussion of the functionality see the paper by Rein & Tamayo 
    2017.

    Requirements
    ------------
    The SimulationArchive Class only works when the following 
    requirements are satisfied.  

    - Only the WHFast, WHFastHelio and IAS15 integrators are supported.
    - Symplectic correcters are supported.
    - The number of particles can not change.
    - Any additional forces or post-timestep modifications that are
      present during the integration also need to be present when 
      the SimulationArchive class is accessed at a later time.


    Examples
    --------
    Here is a simple example:

    >>> sa = rebound.SimulationArchive("archive.bin")
    >>> sim = sa.getSimulation(t=1e6)
    >>> print(sim.particles[1])
    >>> for sim in sa:
    >>>     print(sim.t, sim.particles[1].e)

    """
    _fields_ = [("_inf", c_void_p),
                ("_filename", c_char_p),
                ("version", c_int), 
                ("size_first", c_long), 
                ("size_snapshot", c_long), 
                ("interval", c_double), 
                ("walltime", c_double), 
                ("t0", c_double), 
                ("nblobs", c_long), 
                ("offset", POINTER(c_uint32)), 
                ("t", POINTER(c_double)) 
                ]
    def __init__(self,filename,setup=None, setup_args=(), rebxfilename=None):
        self.setup = setup
        self.setup_args = setup_args
        self.rebxfilename = rebxfilename
        w = c_int(0)
        clibrebound.reb_read_simulationarchive_with_messages(byref(self),c_char_p(filename.encode("ascii")),byref(w))
        if w.value & (1+16+32+64+256) :     # Major error
            raise ValueError(BINARY_WARNINGS[0])
        for message, value in BINARY_WARNINGS:  # Just warnings
            if w.value & value:
                warnings.warn(message, RuntimeWarning)

    """
    def __init__(self,filename,setup=None, setup_args=(), rebxfilename=None):
        self.cfilename = c_char_p(filename.encode("ascii"))
        self.setup = setup
        self.setup_args = setup_args
        self.rebxfilename = rebxfilename

        # Recreate simulation at t=0
        w = c_int(0)
        clibrebound.reb_create_simulation.restype = POINTER(Simulation)
        simp = clibrebound.reb_create_simulation() 
        clibrebound.reb_create_simulation_from_binary_with_messages(simp, self.cfilename,byref(w))
        if (not simp) or (w.value & 1):     # Major error
            raise ValueError(BINARY_WARNINGS[0][0])
        if (w.value & 2):     
            warnings.warn(BINARY_WARNINGS[1][0], RuntimeWarning)
            # Note: Other warnings not shown!
        self.simp = simp
        sim = self.simp.contents
        if self.setup:
            self.setup(sim, *self.setup_args)
        if self.rebxfilename:
            import reboundx
            rebx = reboundx.Extras.from_file(sim, self.rebxfilename)

        self.filesize = os.path.getsize(filename)
        self.dt = sim.dt
        if sim.simulationarchive_interval>0. and sim.simulationarchive_interval_walltime>0.:
            raise ValueError("Something is wrong. Binary file contains both interval>0 and interval_walltime>0.")
        if sim.simulationarchive_interval==0. and sim.simulationarchive_interval_walltime==0.:
            raise ValueError("Something is wrong. Binary file has interval==0 and interval_walltime==0.")
        if sim.simulationarchive_interval>0.:
            self.interval = sim.simulationarchive_interval
        if sim.simulationarchive_interval_walltime>0.:
            self.interval_walltime = sim.simulationarchive_interval_walltime

        self.tmin = sim.t
        clibrebound.reb_simulationarchive_nblobs.restype = c_long
        self.Nblob = clibrebound.reb_simulationarchive_nblobs(simp, self.cfilename)
        if sim.simulationarchive_interval_walltime>0.:
            self.timetable = [-1.]*(self.Nblob+1)
            self.timetable[0] = sim.t
            self._loadAndSynchronize(-1)
            self.timetable[-1] = sim.t
            self.tmax = sim.t
        else:
            self.tmax = self.tmin + self.interval*(self.Nblob)
    """
    
    def __str__(self):
        """
        Returns a string with details of this simulation archive.
        """
        return "<rebound.SimulationArchive instance, snapshots={0} >".format(str(len(self)))

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
        if key>= len(self) or key<0:
            raise IndexError("Index out of range, number of snapshots stored in binary: %d."%len(self))
        
        w = c_int(0)
        sim = Simulation()
        clibrebound.reb_create_simulation_from_simulationarchive_with_messages(byref(sim), byref(self), c_long(key), byref(w))
        if w.value & (1+16+32+64+256) :     # Major error
            raise ValueError(BINARY_WARNINGS[0])
        for message, value in BINARY_WARNINGS:  # Just warnings
            if w.value & value and value!=1:
                warnings.warn(message, RuntimeWarning)
        return sim
    
    def __setitem__(self, key, value):
        raise AttributeError("Cannot modify SimulationArchive.")

    def __delitem__(self, key):
        raise AttributeError("Cannot modify SimulationArchive.")

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        return self.nblobs  # number of SA snapshots (also counting binary at t=0)

    def _loadAndSynchronize(self, snapshot, keep_unsynchronized=1):
        """
        Update self.sim by loading a snapshopt from the binary file. 
        """
        sim = self.simp.contents
        clibrebound.reb_simulationarchive_load_snapshot.restype = c_int
        retv = clibrebound.reb_simulationarchive_load_snapshot(self.simp, self.cfilename, c_long(snapshot))
        if retv:
            raise ValueError("Error while loading snapshot in binary file. Errorcode: %d."%retv)
        if sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1:
            keep_unsynchronized = 0
        if sim.integrator=="mercurius" and sim.ri_mercurius.safe_mode == 1:
            keep_unsynchronized = 0
        sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
        sim.ri_mercurius.keep_unsynchronized = keep_unsynchronized
        sim.integrator_synchronize()
        if snapshot == 0:
            if self.setup:
                self.setup(sim, *self.setup_args)
            if self.rebxfilename:
                import reboundx
                rebx = reboundx.Extras.from_file(sim, self.rebxfilename)

    def _getSnapshotIndex(self, t):
        """
        Return the index for the snapshot just before t
        """
        if t>self.tmax+self.dt or t<self.tmin:
            raise ValueError("Requested time outside of baseline stored in binary fie.")
        try:
            bi = max(0,int(math.floor((t-self.dt-self.tmin)/self.interval)))
            bt = self.tmin + self.interval*bi
            return bi, bt
        except AttributeError: # No interval, need to use timetable
            # Bisection method
            sim = self.simp.contents
            l = 0
            r = self.Nblob
            while True:
                bi = l+(r-l)//2
                if self.timetable[bi] == -1.:
                    clibrebound.reb_simulationarchive_load_snapshot.restype = c_int
                    retv = clibrebound.reb_simulationarchive_load_snapshot(self.simp, self.cfilename, bi)
                    if retv:
                        raise ValueError("Error while loading snapshot in binary file. Errorcode: %d."%retv)
                    self.timetable[bi] = sim.t
                if self.timetable[bi]>t:
                    r = bi
                else:
                    l = bi
                if r-1<=l:
                    bi = l
                    break
            return bi, sim.t
        return 0, 0.

    def getSimulation(self, t, mode='snapshot', keep_unsynchronized=1):
        """
        This function returns a simulation object at (or close to) the requested time `t`. 


        Arguments
        ---------
        t : float
            Requested time. Needs to be within tmin and tmax of this Simulation Archive.
        mode : str
            This argument determines how close the simulation should be to the requested time.
            There are three options. 
            - 'snapshot' This loads a nearby snapshot such that sim.t<t. This is the default.
            - 'close' This integrates the simulation to get to the time t but may overshoot by at most one timestep sim.dt.
            - 'exact' This integrates the simulation to exactly time t. This is not compatible with keep_unsynchronized=1. 
        keep_unsynchronized : int
            By default this argument is 1. This means that if the simulation had to be synchronized to generate this output, then it will nevertheless use the unsynchronized values if one integrates the simulation further in time. This is important for exact (bit-by-bit) reproducibility. If the value of this argument is 0, then one can modify the particles coordinates and these changes are taken into account when integrating the simulation further in time.
        
        Returns
        ------- 
        A rebound.Simulation object. This object will be modified 
        the next time getSimulation is called. Making any manual 
        changes to this object could have unforseen consequences.
        
        Examples
        --------
        Here is a simple example on how to load a simulation from a 
        Simulation Archive file with the `getSimulation` method.
        As the `mode` argument is set to `close`, the simulation
        will be integrated from the nearest snapshot to the request time.

        >>> sa = rebound.SimulationArchive("archive.bin")
        >>> sim = sa.getSimulation(t=1e6, mode="close")
        >>> print(sim.particles[1])
        >>> for sim in sa:
        >>>     print(sim.t, sim.particles[1].e)

        """
        if mode not in ['snapshot', 'close', 'exact']:
            raise AttributeError("Unknown mode.")

        bi, bt = self._getSnapshotIndex(t)
        if mode=='snapshot':
            self._loadAndSynchronize(bi, keep_unsynchronized=keep_unsynchronized)
            return self.simp.contents
        else:
            sim = self.simp.contents
            if sim.t<t and bt-sim.dt<sim.t \
                and ((sim.integrator != "whfast" or (sim.ri_whfast.keep_unsynchronized==1 or sim.ri_whfast.safe_mode == 1))
                or (sim.integrator != "mercurius" or (sim.ri_mercurius.keep_unsynchronized==1 or sim.ri_mercurius.safe_mode == 1))):
                # Reuse current simulation
                pass
            else:
                # Load from snapshot
                clibrebound.reb_simulationarchive_load_snapshot.restype = c_int
                retv = clibrebound.reb_simulationarchive_load_snapshot(self.simp, self.cfilename, c_long(bi));
                if retv:
                    raise ValueError("Error while loading snapshot in binary file. Errorcode: %d."%retv)

            if mode=='exact':
                keep_unsynchronized==0
            if sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1:
                keep_unsynchronized = 0

            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            sim.ri_mercurius.keep_unsynchronized = keep_unsynchronized
            if bi == 0:
                if self.setup:
                    self.setup(sim, *self.setup_args)
                if self.rebxfilename:
                    import reboundx
                    rebx = reboundx.Extras.from_file(sim, self.rebxfilename)

            exact_finish_time = 1 if mode=='exact' else 0
            sim.integrate(t,exact_finish_time=exact_finish_time)
                
            return sim


    def getSimulations(self, times, **kwargs):
        """
        A generator to quickly access many simulations. 
        The arguments are the same as for `getSimulation`.
        """
        for t in times:
            yield self.getSimulation(t, **kwargs)

    
