from ctypes import Structure, c_double, POINTER, c_int, c_uint, c_long, c_ulong, c_void_p, c_char_p, CFUNCTYPE, byref
from . import clibrebound, Escape, NoParticles, Encounter, SimulationError
from .particle import Particle
from .units import units_convert_particle, check_units, convert_G
import math
import os
import ctypes.util
try:
    import pkg_resources
except: 
    # Fails on python3, but not important
    pass
import types
      
### The following enum and class definitions need to
### consitent with those in rebound.h
        
INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "wh": 3, "leapfrog": 4, "hybrid": 5, "none": 6}
BOUNDARIES = {"none": 0, "open": 1, "periodic": 2, "shear": 3}
GRAVITIES = {"none": 0, "basic": 1, "compensated": 2, "tree": 3}
COLLISIONS = {"none": 0, "direct": 1, "tree": 2}

class reb_vec3d(Structure):
    _fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double)]

class reb_dp7(Structure):
    _fields_ = [("p0", POINTER(c_double)),
                ("p1", POINTER(c_double)),
                ("p2", POINTER(c_double)),
                ("p3", POINTER(c_double)),
                ("p4", POINTER(c_double)),
                ("p5", POINTER(c_double)),
                ("p6", POINTER(c_double))]

class reb_ghostbox(Structure):
    _fields_ = [("shiftx", c_double),
                ("shifty", c_double),
                ("shiftz", c_double),
                ("shiftvx", c_double),
                ("shiftvy", c_double),
                ("shiftvz", c_double)]


class reb_simulation_integrator_hybrid(Structure):
    _fields_ = [("switch_ratio", c_double),
                ("mode", c_int)]

class reb_simulation_integrator_wh(Structure):
    _fields_ = [(("allocatedN"), c_int),
                ("eta", POINTER(c_double))]

class reb_simulation_integrator_sei(Structure):
    _fields_ = [("OMEGA", c_double),
                ("OMEGAZ", c_double),
                ("lastdt", c_double),
                ("sindt", c_double),
                ("tandt", c_double),
                ("sindtz", c_double),
                ("tandtz", c_double)]

class reb_simulation_integrator_ias15(Structure):
    _fields_ = [("epsilon", c_double),
                ("min_dt", c_double),
                ("epsilon_global", c_uint),
                ("iterations_max_exceeded", c_ulong),
                ("allocatedN", c_int),
                ("at", POINTER(c_double)),
                ("x0", POINTER(c_double)),
                ("v0", POINTER(c_double)),
                ("a0", POINTER(c_double)),
                ("csx", POINTER(c_double)),
                ("csv", POINTER(c_double)),
                ("csa0", POINTER(c_double)),
                ("g", reb_dp7),
                ("b", reb_dp7),
                ("csb", reb_dp7),
                ("e", reb_dp7),
                ("br", reb_dp7),
                ("er", reb_dp7)]

class reb_simulation_integrator_whfast(Structure):
    _fields_ = [("corrector", c_uint),
                ("recalculate_jacobi_this_timestep", c_uint),
                ("safe_mode", c_uint),
                ("p_j", POINTER(Particle)),
                ("eta", POINTER(c_double)),
                ("Mtotal", c_double),
                ("is_synchronized", c_uint),
                ("allocatedN", c_uint),
                ("timestep_warning", c_uint),
                ("recalculate_jacobi_but_not_synchronized_warning", c_uint)]

class Orbit(Structure):
    """
    A class containing orbital parameters for a particle.
    This is an abstraction of the reb_orbit data structure in C.

    When using the various REBOUND functions using Orbits, all angles are in radians. 
    The following image illustrated the most important angles used.
    In REBOUND the reference direction is the positive x direction, the reference plane
    is the xy plane.
    
    .. image:: images/orbit.png
       :width: 500px
       :height: 450px

    Image from wikipedia. CC-BY-SA-3.

    Attributes
    ----------
    r       : float           
        radial distance from reference 
    v       : float         
        velocity relative to central object's velocity
    h       : float           
        specific angular momentum
    P       : float           
        orbital period
    n       : float           
        mean motion
    a       : float           
        semimajor axis
    e       : float           
        eccentricity
    inc     : float           
        inclination
    Omega   : float           
        longitude of ascending node
    omega   : float           
        argument of pericenter
    pomega  : float           
        longitude of pericenter
    f       : float           
        true anomaly
    M       : float           
        mean anomaly
    l       : float           
        mean longitude = Omega + omega + M
    theta   : float           
        true longitude = Omega + omega + f
    """
    _fields_ = [("r", c_double),
                ("v", c_double),
                ("h", c_double),
                ("P", c_double),
                ("n", c_double),
                ("a", c_double),
                ("e", c_double),
                ("inc", c_double),
                ("Omega", c_double),
                ("omega", c_double),
                ("pomega", c_double),
                ("f", c_double),
                ("M", c_double),
                ("l", c_double),
                ("theta", c_double)]

    def __str__(self):
        """
        Returns a string with the semi-major axis and eccentricity of the orbit.
        """
        return "<rebound.Orbit instance, a={0} e={1} inc={2} Omega={3} omega={4} f={5}>".format(str(self.a),str(self.e), str(self.inc), str(self.Omega), str(self.omega), str(self.f))

class Simulation(Structure):
    """
    REBOUND Simulation Object.

    This object encapsulated an entire REBOUND simulation. 
    It is an abstraction of the C struct reb_simulation.
    You can create mutiple REBOUND simulations (the c library is thread safe). 

    Examples
    --------
    Most simulation parameters can be directly changed with the property syntax:

    >>> sim = rebound.Simulation()
    >>> sim.G = 1.   # Sets the graviational constant (default 1)
    >>> sim.softening = 1.   # Sets the graviational softening parameter default (0)
    >>> sim.dt = 0.1   # Sets the timestep (will change for adaptive integrators such as IAS15).
    >>> sim.t = 0.   # Sets the current simulation time (default 0.)
    >>> print(sim.N)  # Gets the current number of particles
    >>> print(sim.N_active)  # Gets the current number of active particles

    """
    def __init__(self):
        clibrebound.reb_init_simulation(byref(self))

    @classmethod
    def from_file(cls, filename):
        """
        Loads a REBOUND simulation from a file.
        
        After loading the REBOUND simulation from file, you need to reset any function pointers manually.
        
        Arguments
        ---------
        filename : str
            Filename of the binary file.
        
        Returns
        ------- 
        A rebound.Simulation object.
        
        """
        if os.path.isfile(filename):
            clibrebound.reb_create_simulation_from_binary.restype = POINTER_REB_SIM
            return clibrebound.reb_create_simulation_from_binary(c_char_p(filename.encode("ascii"))).contents
        else:
            raise ValueError("File does not exist.")

    def __del__(self):
        if self._b_needsfree_ == 1: # to avoid, e.g., sim.particles[1]._sim.contents.G creating a Simulation instance to get G, and then freeing the C simulation when it immediately goes out of scope
            clibrebound.reb_free_pointers(byref(self))

# Status functions
    def status(self):
        """ 
        Prints a summary of the current status 
        of the simulation.
        """
        from rebound import __version__, __build__
        s= ""
        s += "---------------------------------\n"
        s += "REBOUND version:     \t%s\n" %__version__
        s += "REBOUND built on:    \t%s\n" %__build__
        s += "Number of particles: \t%d\n" %self.N       
        s += "Selected integrator: \t" + self.integrator + "\n"       
        s += "Simulation time:     \t%f\n" %self.t
        s += "Current timestep:    \t%f\n" %self.dt
        if self.N>0:
            s += "---------------------------------\n"
            for p in self.particles:
                s += str(p) + "\n"
        s += "---------------------------------"
        print(s)

# Set function pointer for additional forces
    @property
    def additional_forces(self):
        """
        Get or set a function pointer for calculating additional forces in the simulation.

        The argument can be a python function or something that can 
        be cast to a C function of type CFUNCTYPE(None,POINTER(Simulaton)). 
        If the forces are velocity dependent, the flag 
        force_is_velocity_dependent needs to be set to 1. Otherwise
        the particle structures might contain incorrect velocity 
        values.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @additional_forces.setter
    def additional_forces(self, func):
        self._afp = AFF(func)
        self._additional_forces = self._afp

    @property
    def post_timestep_modifications(self):
        """
        Get or set a function pointer for post-timestep modifications.

        The argument can be a python function or something that can be cast to a C function or a
        python function.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @post_timestep_modifications.setter
    def post_timestep_modifications(self, func):
        self._ptmp = AFF(func)
        self._post_timestep_modifications = self._ptmp
   
    @property 
    def coefficient_of_restitution(self):
        """
        Get or set a function pointer that defined the coefficient of restitution.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @coefficient_of_restitution.setter
    def coefficient_of_restitution(self, func):
        self._corfp = CORFF(func)
        self._coefficient_of_restitution = self._corfp

# Setter/getter of parameters and constants
    @property 
    def N_real(self):
        """
        Get the current number of real (i.e. no variational/shadow) particles in the simulation.
        """
        return self.N-self.N_var

    @property
    def integrator(self):
        """
        Get or set the intergrator module.

        Available integrators are:

        - ``'ias15'`` (default)
        - ``'whfast'``
        - ``'sei'``
        - ``'wh'``
        - ``'leapfrog'``
        - ``'hybrid'``
        - ``'none'``
        
        Check the online documentation for a full description of each of the integrators. 
        """
        i = self._integrator
        for name, _i in INTEGRATORS.items():
            if i==_i:
                return name
        return i
    @integrator.setter
    def integrator(self, value):
        if isinstance(value, int):
            self._integrator = c_int(value)
        elif isinstance(value, basestring):
            debug.integrator_fullname = value
            debug.integrator_package = "REBOUND"
            value = value.lower()
            if value in INTEGRATORS: 
                self._integrator = INTEGRATORS[value]
            elif value.lower() == "mercury":
                debug.integrator_package = "MERCURY"
            elif value.lower() == "swifter-whm":
                debug.integrator_package = "SWIFTER"
            elif value.lower() == "swifter-symba":
                debug.integrator_package = "SWIFTER"
            elif value.lower() == "swifter-helio":
                debug.integrator_package = "SWIFTER"
            elif value.lower() == "swifter-tu4":
                debug.integrator_package = "SWIFTER"
            else:
                raise ValueError("Warning. Integrator not found.")
    
    @property
    def boundary(self):
        """
        Get or set the boundary module.

        Available boundary modules are:

        - ``'none'`` (default)
        - ``'open'``
        - ``'periodic'``
        - ``'shear'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._boundary
        for name, _i in BOUNDARIES.items():
            if i==_i:
                return name
        return i
    @boundary.setter
    def boundary(self, value):
        if isinstance(value, int):
            self._boundary = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in BOUNDARIES: 
                self._boundary = BOUNDARIES[value]
            else:
                raise ValueError("Warning. Boundary condition module not found.")

    @property
    def gravity(self):
        """
        Get or set the gravity module.

        Available gravity modules are:

        - ``'none'`` 
        - ``'basic'``
        - ``'compensated'`` (default)
        - ``'tree'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._gravity
        for name, _i in GRAVITIES.items():
            if i==_i:
                return name
        return i
    @gravity.setter
    def gravity(self, value):
        if isinstance(value, int):
            self._gravity = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in GRAVITIES: 
                self._gravity = GRAVITIES[value]
            else:
                raise ValueError("Warning. Gravity module not found.")

    @property
    def collision(self):
        """
        Get or set the collision module.

        Available collision modules are:

        - ``'none'`` (default)
        - ``'direct'``
        - ``'tree'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._collision
        for name, _i in COLLISIONS.items():
            if i==_i:
                return name
        return i
    @collision.setter
    def collision(self, value):
        if isinstance(value, int):
            self._collision = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in COLLISIONS: 
                self.collision = COLLISIONS[value]
            else:
                raise ValueError("Warning. Collision module not found.")

# Units

    @property
    def units(self):
        """
        Tuple of the units for length, time and mass.  Can be set in any order, and strings are not case-sensitive.  See ipython_examples/Units.ipynb for more information.  You can check the units' exact values and add Additional units in rebound/rebound/units.py.  Units should be set before adding particles to the simulation (will give error otherwise).

        Currently supported Units 
        -------------------------

        Times:
        Hr          : Hours
        Yr          : Julian years
        Jyr         : Julian years
        Sidereal_yr : Sidereal year
        Yr2pi       : Year divided by 2pi, with year defined as orbital period of planet at 1AU around 1Msun star
        Kyr         : Kiloyears (Julian)
        Myr         : Megayears (Julian)
        Gyr         : Gigayears (Julian)

        Lengths:
        M           : Meters
        Cm          : Centimeters
        Km          : Kilometers
        AU          : Astronomical Units

        Masses:
        Kg          : Kilograms
        Msun        : Solar masses
        Mmercury    : Mercury masses
        Mvenus      : Venus masses
        Mearth      : Earth masses
        Mmars       : Mars masses
        Mjupiter    : Jupiter masses
        Msaturn     : Saturn masses
        Muranus     : Neptune masses
        Mpluto      : Pluto masses
        
        Examples
        --------
        
        >>> sim = rebound.Simulation()
        >>> sim.units = ('yr', 'AU', 'Msun')

        """
        if not hasattr(self, '_units'):
            self._units = {'length':None, 'mass':None, 'time':None}

        return self._units

    @units.setter
    def units(self, newunits):
        if not hasattr(self, '_units'):
            self._units = {'length':None, 'mass':None, 'time':None}
        newunits = check_units(newunits)        
        if self.particles: # some particles are loaded
            raise AttributeError("Error:  You cannot set the units after populating the particles array.  See ipython_examples/Units.ipynb.")
        self.update_units(newunits) 

    def update_units(self, newunits): 
        self._units['length'] = newunits[0] 
        self._units['time'] = newunits[1] 
        self._units['mass'] = newunits[2] 
        self.G = convert_G(self._units['length'], self._units['time'], self._units['mass'])

    def convert_particle_units(self, *args): 
        """
        Will convert the units for the simulation (i.e. convert G) as well as the particles' cartesian elements.
        Must have set sim.units ahead of calling this function so REBOUND knows what units to convert from.

        Parameters
        ----------
        3 strings corresponding to units of time, length and mass.  Can be in any order and aren't case sensitive.        You can add new units to rebound/rebound/units.py
        """
        if None in self.units.values():
            raise AttributeError("Must set sim.units before calling convert_particle_units in order to know what units to convert from.")
        new_l, new_t, new_m = check_units(args)
        for p in self.particles:
            units_convert_particle(p, self._units['length'], self._units['time'], self._units['mass'], new_l, new_t, new_m)
        self.update_units((new_l, new_t, new_m))

# MEGNO
    def init_megno(self, delta):
        """
        This function initialises the chaos indicator MEGNO particles and enables their integration.

        MEGNO is short for Mean Exponential Growth of Nearby orbits. It can be used to test
        if a system is chaotic or not. In the backend, the integrator is integrating an additional set
        of particles using the variational equation. Note that variational equations are better 
        suited for this than shadow particles. MEGNO is currently only supported in the IAS15 
        and WHFast integrators.

        This function also needs to be called if you are interested in the Lyapunov exponent as it is
        calculate with the help of MEGNO. See Rein and Tamayo 2015 for details on the implementation.

        The initial delta value can in principle tace any value, typically we choose 1e-16. For 
        more information on MENGO see e.g. http://dx.doi.org/10.1051/0004-6361:20011189
        """
        clibrebound.reb_tools_megno_init(byref(self), c_double(delta))
    
    def calculate_megno(self):
        """
        Return the current MEGNO value.
        Note that you need to call init_megno() before the start of the simulation.
        """
        clibrebound.reb_tools_calculate_megno.restype = c_double
        return clibrebound.reb_tools_calculate_megno(byref(self))
    
    def calculate_lyapunov(self):
        """
        Return the current Lyapunov Characteristic Number (LCN).
        Note that you need to call init_megno() before the start of the simulation.
        To get a timescale (the Lyapunov timescale), take the inverse of this quantity.
        """
        clibrebound.reb_tools_calculate_lyapunov.restype = c_double
        return clibrebound.reb_tools_calculate_lyapunov(byref(self))
    
# Particle add function, used to be called particle_add() and add_particle() 
    def add(self, particle=None, **kwargs):   
        """
        Adds a particle to REBOUND. Accepts one of the following:

        1) A single Particle structure.
        2) The particle's mass and a set of cartesian coordinates: m,x,y,z,vx,vy,vz.
        3) The primary as a Particle structure, the particle's mass and a set of orbital elements primary,a,anom,e,omega,inv,Omega,MEAN (see kepler_particle() for the definition of orbital elements). 
        4) A name of an object (uses NASA Horizons to look up coordinates)
        5) A list of particles or names.
        """
        if particle is not None:
            if isinstance(particle, Particle):
                if kwargs == {}: # copy particle
                    if (self.gravity == "tree" or self.collision == "tree") and self.root_size <=0.:
                        raise ValueError("The tree code for gravity and/or collision detection has been selected. However, the simulation box has not been configured yet. You cannot add particles until the the simulation box has a finite size.")

                    clibrebound.reb_add(byref(self), particle)
                else: # use particle as primary
                    self.add(Particle(simulation=self, primary=particle, **kwargs))
            elif isinstance(particle, list):
                for p in particle:
                    self.add(p)
            elif isinstance(particle,str):
                if None in self.units.values():
                    self.units = ('AU', 'yr2pi', 'Msun')
                self.add(horizons.getParticle(particle,**kwargs))
                units_convert_particle(self.particles[-1], 'km', 's', 'kg', self._units['length'], self._units['time'], self._units['mass'])
        else: 
            self.add(Particle(simulation=self, **kwargs))

# Particle getter functions
    @property
    def particles(self):
        """
        Return an array that points to the particle structure.

        This is an array of pointers and thus the contents of the array update 
        as the simulation progresses. Note that the pointers could change,
        for example when a particle is added or removed from the simulation. 
        """
        ps = []
        N = self.N 
        ps_a = self._particles
        for i in range(0,N):
            ps.append(ps_a[i])
        return ps

    @particles.deleter
    def particles(self):
        """
        Remove all particles from the simulation
        """
        clibrebound.reb_particles_remove_all(byref(self))

    def remove(self, index=None, id=None, keepSorted=1):
        """ 
        Removes a particle from the simulation.

        Parameters
        ----------
        Either the index in the particles array to remove, or the id of the particle to
        remove.  The keepSorted flag ensures the particles array remains sorted
        in order of increasing ids.  One might set it to zero in cases with many
        particles and many removals to speed things up.
        """
        if index is not None:
            success = clibrebound.reb_remove(byref(self), c_int(index), keepSorted)
            if not success:
                print("Index %d passed to remove_particle was out of range (N=%d). Did not remove particle.\n"%(index, self.N))
            return
        if id is not None:
            success = clibrebound.reb_remove_by_id(byref(self), c_int(id), keepSorted)
            if not success:
                print("id %d passed to remove_particle was not found.  Did not remove particle.\n"%(id))


# Orbit calculation
    def calculate_orbits(self, heliocentric=False, barycentric=False):
        """ 
        Calculate orbital parameters for all partices in the simulation.
        By default this functions returns the orbits in Jacobi coordinates. 

        If MEGNO is enabled, variational particles will be ignored.

        Parameters
        ----------
        heliocentric : bool, optional
            Set the parameter heliocentric to True to return orbits referenced to sim.particles[0].
        barycentric : bool, optional
            Set the parameter barycentric to True to return orbits referenced to the system's barycenter.
            

        Returns
        -------
        Returns an array of Orbits of length N-1.
        """
        _particles_tmp = self.particles
        orbits = []
        
        jacobi = True
        com = _particles_tmp[0]
        if heliocentric is True:
            jacobi = False
        if barycentric is True:
            com = self.calculate_com()
            jacobi = False

        clibrebound.reb_get_com_of_pair.restype = Particle
        for i in range(1,self.N_real):
            orbits.append(_particles_tmp[i].calculate_orbit(primary=com))
            if jacobi is True:
                com = clibrebound.reb_get_com_of_pair(com, _particles_tmp[i])

        return orbits

# COM calculation 
    def calculate_com(self, last=None):
        """
        Returns the center of momentum for all particles in the simulation.

        Parameters
        ----------
        last : int or None, optional
            If ``last`` is specified only calculate the center of momentum for the
            first ``last`` particles in the array (i.e., indices up to i-1, as used 
            in Jacobi coordinates).

        Examples
        --------
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1, x=0)
        >>> sim.add(m=1, x=1)
        >>> com = sim.calculate_com()
        >>> com.x
        0.5

        """
        if last is not None:
            last = min(last, self.N_real-1)
            clibrebound.reb_get_jacobi_com.restype = Particle
            com = clibrebound.reb_get_jacobi_com(byref(self.particles[last]))
            return com
        else:
            clibrebound.reb_get_com.restype = Particle
            com = clibrebound.reb_get_com(byref(self))
            return com
        

# Tools
    def move_to_com(self):
        """
        This function moves all particles in the simulation to a center of momentum frame.
        In that frame, the center of mass is at the origin and does not move.
        It makes sense to call this function at the beginning of the integration, especially 
        for the high accuray integrators IAS15 and WHFast.
        """
        clibrebound.reb_move_to_com(byref(self))
    
    def calculate_energy(self):
        """
        Returns the sum of potential and kinetic energy of all particles in the simulation.
        """
        clibrebound.reb_tools_energy.restype = c_double
        return clibrebound.reb_tools_energy(byref(self))
    
    def configure_box(self, boxsize, root_nx=1, root_ny=1, root_nz=1):
        """
        Initialize the simulation box.

        This function only needs to be called it boundary conditions other than "none"
        are used. In such a case the boxsize must be known and is set with this function.

        Parameters
        ----------
        boxsize : float, optional
            The size of one root box.
        root_nx, root_ny, root_nz : int, optional
            The number of root boxes in each direction. The total size of the simulation box
            will be ``root_nx * boxsize``, ``root_ny * boxsize`` and ``root_nz * boxsize``.
            By default there will be exactly one root box in each direction.
        """
        clibrebound.reb_configure_box(byref(self), c_double(boxsize), c_int(root_nx), c_int(root_ny), c_int(root_nz))
        return
   
    def configure_ghostboxes(self, nghostx=0, nghosty=0, nghostz=0):
        """
        Initialize the ghost boxes.

        This function only needs to be called it boundary conditions other than "none" or
        "open" are used. In such a case the number of ghostboxes must be known and is set 
        with this function. 
        
        Parameters
        ----------
        nghostx, nghosty, nghostz : int
            The number of ghost boxes in each direction. All values default to 0 (no ghost boxes).
        """
        clibrebound.nghostx = c_int(nghostx)
        clibrebound.nghosty = c_int(nghosty)
        clibrebound.nghostz = c_int(nghostz)
        return


# Input/Output routines
    def save(self, filename):
        """
        Save the entire REBOUND simulation to a binary file.
        """
        clibrebound.reb_output_binary(byref(self), c_char_p(filename.encode("ascii")))
        
# Integrator Flags
    @property 
    def integrator_sei_OMEGA(self):
        """
        Get or set the epicyclic frequency OMEGA.

        For simulations making use of shearing sheet boundary conditions,
        REBOUND needs to know the epicyclic frequency. By default OMEGA
        is 1. For more details read Rein and Tremaine 2011.
        """
        return self.ri_sei.OMEGA
    @integrator_sei_OMEGA.setter 
    def integrator_sei_OMEGA(self, value):
        self.ri_sei.OMEGA = c_double(value)

    @property 
    def integrator_whfast_corrector(self):
        """
        Get or set the order of the symplectic corrector in the WHFast integrator.

        By default the symplectic correctors are turned off (=0). For high
        accuracy simulation set this value to 11. For more details read 
        Rein and Tamayo 2015.
        """
        return self.ri_whfast.corrector
    @integrator_whfast_corrector.setter 
    def integrator_whfast_corrector(self, value):
        self.ri_whfast.corrector = c_uint(value)

    @property
    def integrator_whfast_safe_mode(self):
        """
        Get or set the safe mode flag for WHFast.

        If safe_mode is 1 (default) particles can be modified between
        timesteps and particle velocities and positions are always synchronised.
        If you set safe_mode to 0, the speed and accuracy of WHFast improves.
        However, make sure you are aware of the consequences. Read the iPython tutorial
        on advanced WHFast usage to learn more.
        """
        return self.ri_whfast.safe_mode
    @integrator_whfast_safe_mode.setter
    def integrator_whfast_safe_mode(self, value):
        self.ri_whfast.safe_mode = c_uint(value)

    @property
    def integrator_whfast_recalculate_jacobi_this_timestep(self):
        """
        Sets a flag that tells WHFast that the particles have changed.

        Setting this flag to 1 (default 0) triggers the WHFast integrator
        to recalculate Jacobi coordinates. This is needed if the user changes 
        the particle position, velocity or mass inbetween timesteps.
        After every timestep the flag is set back to 0, so if you continuously
        update the particles manually, you need to set this flag to 1 after every timestep.
        """ 
        return self.ri_whfast.recalculate_jacobi_this_timestep
    @integrator_whfast_recalculate_jacobi_this_timestep.setter
    def integrator_whfast_recalculate_jacobi_this_timestep(self, value):
        self.ri_whfast.recalculate_jacobi_this_timestep = c_uint(value)
    
# Integration
    def step(self):
        """
        Perform exactly one integration step with REBOUND. This function is rarely needed.
        Instead, use integrate().
        """
        clibrebound.reb_step(byref(self))

    def integrate(self, tmax, exact_finish_time=1):
        """
        Main integration function. Call this function when you have setup your simulation and want to integrate it forward (or backward) in time. The function might be called many times to integrate the simulation in steps and create outputs in-between steps.
        
        Parameters
        ----------
        tmax : float
            The final time of your simulation. If the current time is 100, and tmax=200, then after the calling the integrate routine, the time has advanced to t=200. If tmax is larger than or equal to the current time, no integration will be performed.
        exact_finish_time: int, optional
            This argument determines whether REBOUND should try to finish at the exact time (tmax) you give it or if it is allowed to overshoot. Overshooting could happen if one starts at t=0, has a timestep of dt=10 and wants to integrate to tmax=25. With ``exact_finish_time=1``, the integrator will choose the last timestep such that t is exactly 25 after the integration, otherwise t=30. Note that changing the timestep does affect the accuracy of symplectic integrators negatively.
        
        Exceptions
        ----------
        Exceptions are thrown when no more particles are left in the simulation or when a generic integration error occured. 
        If you specified exit_min_distance or exit_max_distance, then additional exceptions might thrown for escaping particles or particles that undergo a clos encounter.
        
        Examples
        -------- 
        The typical usage is as follows. Note the use of ``np.linspace`` to create equally spaced outputs. 
        Using ``np.logspace`` can be used to easily produce logarithmically spaced outputs.

        >>> import numpy as np
        >>> for time in np.linspace(0,100.,10): 
        >>>     sim.integrate(time)
        >>>     perform_output(sim)
        
        """
        if debug.integrator_package =="REBOUND":
            clibrebound.reb_integrate.restype = c_int
            self.exact_finish_time = c_int(exact_finish_time)
            ret_value = clibrebound.reb_integrate(byref(self), c_double(tmax))
            if ret_value == 1:
                raise SimulationError("An error occured during the integration.")
            if ret_value == 2:
                raise NoParticles("No more particles left in simulation.")
            if ret_value == 3:
                raise Encounter("Two particles had a close encounter (d<exit_min_distance).")
            if ret_value == 4:
                raise Escape("A particle escaped (r>exit_max_distance).")
        else:
            debug.integrate_other_package(tmax,exact_finish_time)

    def integrator_synchronize(self):
        """
        Call this function if safe-mode is disabled and you need synchronize particle positions and velocities between timesteps.
        """
        clibrebound.reb_integrator_synchronize(byref(self))


# Setting up fields after class definition (because of self-reference)
Simulation._fields_ = [
                ("t", c_double),
                ("G", c_double),
                ("softening", c_double),
                ("dt", c_double),
                ("dt_last_done", c_double),
                ("N", c_int),
                ("N_var", c_int),
                ("N_active", c_int),
                ("allocated_N", c_int),
                ("_particles", POINTER(Particle)),
                ("gravity_cs", POINTER(reb_vec3d)),
                ("gravity_cs_allocatedN", c_int),
                ("tree_root", c_void_p),
                ("tree_needs_update", c_int),
                ("opening_angle2", c_double),
                ("_status", c_int),
                ("exact_finish_time", c_int),
                ("force_is_velocity_dependent", c_uint),
                ("gravity_ignore_10", c_uint),
                ("output_timing_last", c_double),
                ("exit_max_distance", c_double),
                ("exit_min_distance", c_double),
                ("usleep", c_double),
                ("boxsize", reb_vec3d),
                ("boxsize_max", c_double),
                ("root_size", c_double),
                ("root_n", c_int),
                ("root_nx", c_int),
                ("root_ny", c_int),
                ("root_nz", c_int),
                ("nghostx", c_int),
                ("nghosty", c_int),
                ("nghostz", c_int),
                ("collisions", c_void_p),
                ("collisions_allocatedN", c_int),
                ("minimum_collision_celocity", c_double),
                ("collisions_plog", c_double),
                ("max_radius", c_double*2),
                ("collisions_Nlog", c_long),
                ("_calculate_megno", c_int),
                ("megno_Ys", c_double),
                ("megno_Yss", c_double),
                ("megno_cov_Yt", c_double),
                ("megno_var_t", c_double),
                ("megno_mean_t", c_double),
                ("megno_mean_Y", c_double),
                ("megno_n", c_long),
                ("_collision", c_int),
                ("_integrator", c_int),
                ("_boundary", c_int),
                ("_gravity", c_int),
                ("ri_sei", reb_simulation_integrator_sei), 
                ("ri_wh", reb_simulation_integrator_wh), 
                ("ri_hybrid", reb_simulation_integrator_hybrid),
                ("ri_whfast", reb_simulation_integrator_whfast),
                ("ri_ias15", reb_simulation_integrator_ias15),
                ("_additional_forces", CFUNCTYPE(None,POINTER(Simulation))),
                ("_post_timestep_modifications", CFUNCTYPE(None,POINTER(Simulation))),
                ("_heartbeat", CFUNCTYPE(None,POINTER(Simulation))),
                ("_coefficient_of_restitution", CFUNCTYPE(c_double,POINTER(Simulation), c_double)),
                ("_collisions_resolve", CFUNCTYPE(None,POINTER(Simulation), c_void_p)),
                ("extras", c_void_p),
                 ]

Particle._fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double),
                ("vx", c_double),
                ("vy", c_double),
                ("vz", c_double),
                ("ax", c_double),
                ("ay", c_double),
                ("az", c_double),
                ("m", c_double),
                ("r", c_double),
                ("lastcollision", c_double),
                ("c", c_void_p),
                ("id", c_int),
                ("ap", c_void_p),
                ("_sim", POINTER(Simulation))]

POINTER_REB_SIM = POINTER(Simulation) 
AFF = CFUNCTYPE(None,POINTER_REB_SIM)
CORFF = CFUNCTYPE(c_double,POINTER_REB_SIM, c_double)

# Import at the end to avoid circular dependence
from . import horizons
from . import debug
