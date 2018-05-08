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

    The SimulationArchive is a binary file format which includes all
    settings, constants as well as particle positions and velocities.
    This makes it possible to reproduce a simulation exactly 
    (down to the last bit). The SimulationArchive allows you to add
    an arbitrary number of snapshots. Simulations can be reconstructed
    from these snapshots. Since version 2 of the SimulationArchive
    (Spring 2018), you can change anything inbetween snapshots,
    icluding settings like the integrator, the timestep, the number of
    particles. The file format is efficient in that only data
    that changed is stored in the SimulationArchive file. This is all
    done automatically. All the user has to do is call the function
    to create a snapshot.
    The SimulationArchive thus allows for fast access to any long running 
    simulations. For a full discussion of the functionality see the paper 
    by Rein & Tamayo 2017.

    Requirements
    ------------
    When using the SimulationArchive, the user is responsible for 
    setting up any additional forces or post-timestep modifications that 
    were  present during the original integration. 
        
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
                ("auto_interval", c_double), 
                ("auto_walltime", c_double), 
                ("nblobs", c_long), 
                ("offset", POINTER(c_uint32)), 
                ("t", POINTER(c_double)) 
                ]
    def __init__(self,filename,setup=None, setup_args=(), rebxfilename=None):
        """
        Arguments
        ---------
        filename : str
            Filename of the SimulationArchive file to be opened.
        setup : function
            Function to be called everytime a simulation object is created
            In this function, the user can setup additional forces
        setup_args : list
            Arguments passed to setup function.
        rebxfilename : str
            Filename of the REBOUNDx binary file.

        """
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
        if self.nblobs<1:
            RuntimeError("Something went wrong. SimulationArchive is empty.")
        self.tmin = self.t[0]
        self.tmax = self.t[self.nblobs-1]

    def __del__(self):
        if self._b_needsfree_ == 1: 
            clibrebound.reb_free_simulationarchive_pointers(byref(self))

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
        if self.setup:
            self.setup(sim, *self.setup_args)
        if self.rebxfilename:
            import reboundx
            rebx = reboundx.Extras.from_file(sim, self.rebxfilename)
        if w.value & (1+16+32+64+256) :     # Major error
            raise ValueError(BINARY_WARNINGS[0])
        for message, value in BINARY_WARNINGS:  # Just warnings
            if w.value & value:
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

    def _getSnapshotIndex(self, t):
        """
        Return the index for the snapshot just before t
        """
        if t>self.tmax or t<self.tmin:
            raise ValueError("Requested time outside of baseline stored in binary file.")
        # Bisection method
        l = 0
        r = len(self)
        while True:
            bi = l+(r-l)//2
            if self.t[bi]>t:
                r = bi
            else:
                l = bi
            if r-1<=l:
                bi = l
                break
        return bi, self.t[bi]

    def getSimulation(self, t, mode='snapshot', keep_unsynchronized=1):
        """
        This function returns a simulation object at (or close to) the requested time `t`. 
        Everytime this function is called a new simulation object is created.


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
        A rebound.Simulation object. Everytime the function gets called
        a new object gets created. 
        
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
        sim = Simulation()
        w = c_int(0)
        clibrebound.reb_create_simulation_from_simulationarchive_with_messages(byref(sim),byref(self),bi,byref(w))

        # Restore function pointers and any additional setup required by the user/reboundx provided functions
        if self.setup:
            self.setup(sim, *self.setup_args)
        if self.rebxfilename:
            import reboundx
            rebx = reboundx.Extras.from_file(sim, self.rebxfilename)

        if mode=='snapshot':
            if sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1:
                keep_unsynchronized = 0
            if sim.integrator=="mercurius" and sim.ri_mercurius.safe_mode == 1:
                keep_unsynchronized = 0
            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            sim.ri_mercurius.keep_unsynchronized = keep_unsynchronized
            sim.integrator_synchronize()
            return sim
        else:
            if mode=='exact':
                keep_unsynchronized==0
            if sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1:
                keep_unsynchronized = 0

            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            sim.ri_mercurius.keep_unsynchronized = keep_unsynchronized
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

    
