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
                ("auto_step", c_ulonglong), 
                ("nblobs", c_long), 
                ("offset", POINTER(c_uint32)), 
                ("t", POINTER(c_double)) 
                ]
    def __init__(self,filename,setup=None, setup_args=(), rebxfilename=None, process_warnings=True):
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
        process_warnings : Bool
            Display warning messages if True (default). Only fail on major errors if set to False.

        """
        self.setup = setup
        self.setup_args = setup_args
        self.rebxfilename = rebxfilename
        w = c_int(0)
        clibrebound.reb_read_simulationarchive_with_messages(byref(self),c_char_p(filename.encode("ascii")),byref(w))
        for majorerror, value, message in BINARY_WARNINGS:
            if w.value & value:
                if majorerror:
                    raise RuntimeError(message)
                else:  
                    # Just a warning
                    if process_warnings:
                        warnings.warn(message, RuntimeWarning)
        else:
            # Store for later
            self.warnings = w
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
        for majorerror, value, message in BINARY_WARNINGS:
            if w.value & value:
                if majorerror:
                    raise RuntimeError(message)
                else:  
                    # Just a warning
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
            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            sim.integrator_synchronize()
            return sim
        else:
            if mode=='exact':
                keep_unsynchronized==0
            if sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1:
                keep_unsynchronized = 0

            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
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

    
    def getBezierPaths(self,origin=None):
        """
        This function returns array that can be used as a Cubic Bezier
        Path in matplotlib. 
        The function returns two arrays, the first one contains
        the verticies for each particles and has the shape
        (Nvert, Nparticles, 2) where Nvert is the number of verticies.
        The second array returned describes the type of verticies to be
        used with matplotlib's Patch class.

        Arguments
        ---------
        origin : multiple, optional
                 If `origin` is None (default), then none of the 
                 coordinates are shifted. If `origin` is an integer
                 then the particle with that index is used as the
                 origin. if `origin` is equal to `com`, then the 
                 centre of mass is used as the origin. 


        Examples
        --------
        The following example reads in a SimulationArchive and plots
        the trajectories as Cubic Bezier Curves. It also plots the 
        actual datapoints stored in the SimulationArchive. 
        Note that the SimulationArchive needs to have enough
        datapoints to allow for smooth and reasonable orbits.

        >>> from matplotlib.path import Path
        >>> import matplotlib.patches as patches
        >>> sa = rebound.SimulationArchive("test.bin")
        >>> verts, codes = sa.getBezierPaths(origin=0)
        >>> fig, ax = plt.subplots()
        >>> for j in range(sa[0].N):
        >>>     path = Path(verts[:,j,:], codes)
        >>>     patch = patches.PathPatch(path, facecolor='none')
        >>>     ax.add_patch(patch)
        >>>     ax.scatter(verts[::3,j,0],verts[::3,j,1])
        >>> ax.set_aspect('equal')
        >>> ax.autoscale_view()
        
        """
        import numpy as np
        Npoints = len(self)*3-2
        if len(self)<=1:
            raise Runtim
        Nparticles = self[0].N
        verts = np.zeros((Npoints,Nparticles,2))
        xy = np.zeros((len(self),Nparticles,2))

        if origin=="com":
            origin = -2
        elif origin is not None:
            try:
                origin = int(origin)
            except:
                raise AttributeError("Cannot parse origin")
            if origin<0 or origin>=Nparticles:
                raise AttributeError("Origin index out of range")


        for i, sim in enumerate(self):
            if origin is None:
                shift = (0,0,0,0)
            elif origin == -2:
                sp = sim.calculate_com()
                shift = (sp.x,sp.y,sp.vx,sp.vy)
            else:
                sp = sim.particles[origin]
                shift = (sp.x,sp.y,sp.vx,sp.vy)
            for j in range(sim.N):
                p = sim.particles[j]
                if i==0:
                    verts[0,j] = p.x-shift[0],p.y-shift[1]
                    verts[1,j] = p.vx-shift[2], p.vy-shift[3]
                else:
                    dt = sim.t-tlast # time since last snapshot
                    verts[-2+i*3,j] = verts[-2+i*3,j]*dt/3.+verts[-3+i*3,j]

                    verts[ 0+i*3,j] = p.x-shift[0],p.y-shift[1]

                    verts[-1+i*3,j] = -p.vx+shift[2], -p.vy+shift[3]
                    verts[-1+i*3,j] = verts[-1+i*3+0,j]*dt/3.+verts[ 0+i*3,j]

                    if i!=len(self)-1:
                        verts[+1+i*3,j] = p.vx-shift[2], p.vy-shift[3]


                xy[i,j] = p.x,p.y
            tlast = sim.t
        codes = np.full(Npoints,4,dtype=np.uint8) # Hardcoded 4 = matplotlib.path.Path.CURVE4
        codes[0] = 1 # Hardcoded 1 = matplotlib.path.Path.MOVETO
        return verts, codes
