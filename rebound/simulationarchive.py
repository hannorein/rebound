from ctypes import Structure, c_double, POINTER, c_int, c_int64, c_uint64, c_void_p, c_char_p, byref 
import os
import sys
import math
import warnings

class Simulationarchive(Structure):
    """
    Simulationarchive Class.

    The Simulationarchive is a binary file format which includes all
    settings, constants as well as particle positions and velocities.
    This makes it possible to reproduce a simulation exactly 
    (down to the last bit). The Simulationarchive allows you to add
    an arbitrary number of snapshots. Simulations can be reconstructed
    from these snapshots. Since version 2 of the Simulationarchive
    (Spring 2018), you can change anything in-between snapshots,
    including settings like the integrator, the timestep, the number of
    particles. The file format is efficient in that only data
    that changed is stored in the Simulationarchive file. This is all
    done automatically. All the user has to do is call the function
    to create a snapshot.
    The Simulationarchive thus allows for fast access to any long-running 
    simulations. For a full discussion of the functionality see the paper 
    by Rein & Tamayo 2017.

    Requirements
    ------------
    When using the Simulationarchive, the user is responsible for 
    setting up any additional forces or post-timestep modifications that 
    were  present during the original integration. 
        
    Examples
    --------
    Here is a simple example:

    >>> sa = rebound.Simulationarchive("archive.bin")
    >>> sim = sa.getSimulation(t=1e6)
    >>> print(sim.particles[1])
    >>> for sim in sa:
    >>>     print(sim.t, sim.particles[1].e)

    """
    _fields_ = [("_inf", c_void_p),
                ("_filename", c_char_p),
                ("version", c_int), 
                ("_reb_version_major", c_int), 
                ("_reb_version_minor", c_int), 
                ("_reb_version_patch", c_int), 
                ("auto_interval", c_double), 
                ("auto_walltime", c_double), 
                ("auto_step", c_uint64), 
                ("nblobs", c_int64), 
                ("offset", POINTER(c_uint64)), 
                ("t", POINTER(c_double)) 
                ]
    def __repr__(self):
        return '<{0}.{1} object at {2}, nblobs={3}, reb_version={4}.{5}.{6}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.nblobs, self._reb_version_major, self._reb_version_minor, self._reb_version_patch)

    def __init__(self,filename,setup=None, setup_args=(), process_warnings=True, reuse_index=None):
        """
        Arguments
        ---------
        filename : str or bytes
            Filename of the Simulationarchive file to be opened.
            Can also be of type bytes to read from memory (uses fmemopen).
        setup : function
            Function to be called everytime a simulation object is created
            In this function, the user can setup additional forces
        setup_args : list
            Arguments passed to setup function.
        process_warnings : Bool
            Display warning messages if True (default). Only fail on major errors if set to False.
        reuse_index : Simulationarchive
            Useful when loading many large Simulationarchives. After loading the first 
            Simulationarchive, pass it as this argument when opening other Simulationarchives with the 
            same shape. Note: Simulationarchive shape must be exactly the same to avoid unexpected
            behaviour.

        """
        self.setup = setup
        self.setup_args = setup_args
        self.process_warnings = process_warnings
        w = c_int(0)
        if reuse_index:
            # Optimized loading
            clibrebound.reb_simulationarchive_create_from_file_with_messages(byref(self),c_char_p(filename.encode("ascii")), byref(reuse_index), byref(w))

        else:
            clibrebound.reb_simulationarchive_create_from_file_with_messages(byref(self),c_char_p(filename.encode("ascii")), None, byref(w))
        for majorerror, value, message in BINARY_WARNINGS:
            if w.value & value:
                if majorerror:
                    raise RuntimeError(message)
                else:  
                    # Just a warning
                    if value==2: # Version warning. Append version used to save SA to message
                        sa_version = "%d.%d.%d" %(self._reb_version_major, self._reb_version_minor, self._reb_version_patch)
                        if sa_version != __version__ and sa_version != "0.0.0":
                            message += " Binary file was saved with REBOUND Version " + sa_version + "."
                            message += " You are currently using REBOUND Version " +  __version__ + "."
                    if process_warnings:
                        warnings.warn(message, RuntimeWarning)
        if not process_warnings:
            # Store for later
            self.warnings = w
        if self.nblobs<1:
            RuntimeError("Something went wrong. Simulationarchive is empty.")
        self.tmin = self.t[0]
        self.tmax = self.t[self.nblobs-1]

    def __del__(self):
        if self._b_needsfree_ == 1: 
            clibrebound.reb_simulationarchive_free_pointers(byref(self))

    def __str__(self):
        """
        Returns a string with details of this simulationarchive.
        """
        return '<{0}.{1} object at {2}, nblobs={3}, reb_version={4}.{5}.{6}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.nblobs, self._reb_version_major, self._reb_version_minor, self._reb_version_patch)

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
        clibrebound.reb_simulation_create_from_simulationarchive_with_messages(byref(sim), byref(self), c_int64(key), byref(w))
        if self.setup:
            self.setup(sim, *self.setup_args)
        for majorerror, value, message in BINARY_WARNINGS:
            if w.value & value:
                if majorerror:
                    raise RuntimeError(message)
                else:  
                    # Just a warning
                    if self.process_warnings:
                        warnings.warn(message, RuntimeWarning)
        if sim.ri_eos.is_synchronized==0 or sim.ri_mercurius.is_synchronized==0 or sim.ri_whfast.is_synchronized==0 or sim.ri_mercurius.is_synchronized==0:
            warnings.warn("The simulation might not be synchronized. You can manually synchronize it by calling sim.synchronize().", RuntimeWarning)

        return sim
    
    def __setitem__(self, key, value):
        raise AttributeError("Cannot modify Simulationarchive.")

    def __delitem__(self, key):
        raise AttributeError("Cannot modify Simulationarchive.")

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
            Requested time. Needs to be within tmin and tmax of this Simulationarchive.
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
        Simulationarchive file with the `getSimulation` method.
        As the `mode` argument is set to `close`, the simulation
        will be integrated from the nearest snapshot to the request time.

        >>> sa = rebound.Simulationarchive("archive.bin")
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
        clibrebound.reb_simulation_create_from_simulationarchive_with_messages(byref(sim),byref(self),bi,byref(w))

        # Restore function pointers and any additional setup required by the user provided functions
        if self.setup:
            self.setup(sim, *self.setup_args)

        if mode=='snapshot':
            if (sim.integrator=="mercurius" and sim.ri_mercurius.safe_mode == 1) or (sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1) or (sim.integrator=="saba" and sim.ri_saba.safe_mode == 1):
                keep_unsynchronized = 0
            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            sim.ri_saba.keep_unsynchronized = keep_unsynchronized
            sim.synchronize()
            return sim
        else:
            if mode=='exact':
                keep_unsynchronized = 0
            if (sim.integrator=="mercurius" and sim.ri_mercurius.safe_mode == 1) or (sim.integrator=="whfast" and sim.ri_whfast.safe_mode == 1) or (sim.integrator=="saba" and sim.ri_saba.safe_mode == 1):
                keep_unsynchronized = 0
            sim.ri_whfast.keep_unsynchronized = keep_unsynchronized
            sim.ri_saba.keep_unsynchronized = keep_unsynchronized
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
        the vertices for each particle and has the shape
        (Nvert, Nparticles, 2) where Nvert is the number of vertices.
        The second array returned describes the type of vertices to be
        used with matplotlib's Patch class.

        Arguments
        ---------
        origin : multiple, optional
                 If `origin` is None (default), then none of the 
                 coordinates are shifted. If `origin` is an integer
                 then the particle with that index is used as the
                 origin. if `origin` is equal to `com`, then the 
                 center of mass is used as the origin. 


        Examples
        --------
        The following example reads in a Simulationarchive and plots
        the trajectories as Cubic Bezier Curves. It also plots the 
        actual datapoints stored in the Simulationarchive. 
        Note that the Simulationarchive needs to have enough
        datapoints to allow for smooth and reasonable orbits.

        >>> from matplotlib.path import Path
        >>> import matplotlib.patches as patches
        >>> sa = rebound.Simulationarchive("test.bin")
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
            raise RuntimeError()
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
                sp = sim.com()
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

from .simulation import Simulation, BINARY_WARNINGS
from . import clibrebound, __version__
