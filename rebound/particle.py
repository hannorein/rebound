from ctypes import Structure, c_double, c_int, byref, memmove, sizeof, c_uint32, c_uint, c_uint64, string_at, POINTER, c_char, c_void_p
import math
import sys

def notNone(a):
    """
    Returns True if array a contains at least one element that is not None. Returns False otherwise.
    """
    return a.count(None) != len(a)

class Particle(Structure):
    """
    The main REBOUND particle data structure. 
    This is an abstraction of the reb_particle structure in C.
    The Particle fields are set at the end of simulation.py to avoid circular references.
    
    Attributes
    ----------
    x, y, z     : float       
        Particle positions
    vx, vy, vz  : float       
        Particle velocities
    ax, ay, az  : float       
        Particle accelerations
    m           : float       
        Particle mass
    r           : float       
        Particle radius
    last_collision : float       
        Last time the particle had a physical collision (if checking for collisions)
    c           : c_void_p (C void pointer) 
        Pointer to the cell the particle is currently in (if using tree code)
    hash          : c_uint32         
        Particle hash (permanent identifier for the particle)
    ap          : c_void_p (C void pointer)
        Pointer to additional parameters one might want to add to particles
    _sim        : POINTER(rebound.Simulation)
        Internal pointer to the parent simulation (used in C version of REBOUND)
    a, e, inc, Omega, omega, f	: float
	    (Kepler Elements) Semi-major axis, eccentricity, inclination, longitude of the ascending node, argument of periapsis, and true anomaly respectively. The Keplerian Elements are in Jacobi coordinates (with mu = G*Minc, where Minc is the total mass from index 0 to the particle's index, inclusive).
    """
    def __repr__(self):
        """ 
        Returns a string with the position and velocity of the particle.
        """
        return '<{0}.{1} object at {2}, m={3} x={4} y={5} z={6} vx={7} vy={8} vz={9}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.m, self.x, self.y, self.z, self.vx, self.vy, self.vz)
   

    def __init__(self, simulation=None, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, P=None, e=None, inc=None, Omega=None, omega=None, pomega=None, f=None, M=None, E=None, l=None, theta=None, T=None, r=None, date=None, variation=None, variation2=None, h=None, k=None, ix=None, iy=None, pal_h=None, pal_k=None, pal_ix=None, pal_iy=None, hash=0, jacobi_masses=False):
        """
        Initializes a Particle structure. Rather than explicitly creating 
        a Particle structure, users may use the ``add()`` member function 
        of a Simulation instance, which will both create a Particle and 
        then add it to the simulation with one function call.

        This function accepts either cartesian positions and velocities, 
        classical orbital elements together with the reference Particle 
        (the primary), as well as orbital parameters defined by Pal (2009).

        For convenience, optional keywords that are not passed default 
        to zero (mass, cartesian and orbital elements). 

        Whenever initializing a particle from orbital elements, one must 
        specify either the semimajor axis or the period of the orbit.
        
        For classical orbital parameters, one can specify the longitude 
        of the ascending node by passing Omega, to specify the pericenter 
        one can pass either omega or pomega (not both), and for the 
        longitude/anomaly one can pass one of f, M, l or theta.  
        See ipython_examples/OrbitalElements.ipynb for examples.  
        See also Murray & Dermott Solar System Dynamics for formal 
        definitions of angles in orbital mechanics.

        All angles should be specified in radians.
        

        Parameters
        ----------
        simulation  : Simulation  
            Simulation instance associated with this particle (Required if passing orbital elements or setting up a variation).
        particle    : Particle, optional    
            If a particle is passed, a copy of that particle is returned.
            If a variational particle is initialized, then ``particle`` is 
            original particle that will be varied. 
        m           : float       
            Mass        (Default: 0)
        x, y, z     : float       
            Positions in Cartesian coordinates  (Default: 0)
        vx, vy, vz  : float       
            Velocities in Cartesian coordinates (Default: 0)
        primary     : Particle    
            Primary body for converting orbital elements to cartesian (Default: center of mass of the particles in the passed simulation, i.e., this will yield Jacobi coordinates as one progressively adds particles) 
        a           : float       
            Semimajor axis (a or P required if passing orbital elements)
        P           : float
            Orbital period (a or P required if passing orbital elements)
        e           : float       
            Eccentricity                (Default: 0)
        inc         : float       
            Inclination                 (Default: 0)
        Omega       : float       
            Longitude of ascending node (Default: 0)
        omega       : float       
            Argument of pericenter      (Default: 0)
        pomega      : float       
            Longitude of pericenter     (Default: 0)
        f           : float       
            True anomaly                (Default: 0)
        M           : float       
            Mean anomaly                (Default: 0)
        E           : float       
            Eccentric anomaly           (Default: 0)
        l           : float       
            Mean longitude              (Default: 0)
        theta       : float       
            True longitude              (Default: 0)
        T           : float 
            Time of pericenter passage  
        h or pal_h  : float       
            h variable, see Pal (2009) for a definition  (Default: 0)
        k or pal_k  : float       
            k variable, see Pal (2009) for a definition  (Default: 0)
        ix or pal_ix : float       
            ix variable, see Pal (2009) for a definition  (Default: 0)
        iy or pal_iy : float       
            iy variable, see Pal (2009) for a definition  (Default: 0)
        r           : float       
            Particle radius (only used for collisional simulations)
        date        : string      
            For consistency with adding particles through horizons.  Not used here.
        variation   : string            (Default: None)
            Set this string to the name of an orbital parameter to initialize the particle as a variational particle.
            Can be one of the following: m, a, e, inc, omega, Omega, f, k, h, lambda, ix, iy.
        variation2  : string            (Default: None)
            Set this string to the name of a second orbital parameter to initialize the particle as a second order variational particle. Only used for second order variational equations. 
            Can be one of the following: m, a, e, inc, omega, Omega, f, k, h, lambda, ix, iy.
        hash        : c_uint32  
            Unsigned integer identifier for particle.  Can pass an integer directly, or a string that will be converted to a hash. User is responsible for assigning unique hashes.
        jacobi_masses: bool
            Whether to use jacobi primary mass in orbit initialization. Particle mass will still be set to physical value (Default: False)
        Examples
        --------

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> p1 = rebound.Particle(simulation=sim, m=0.001, a=0.5, e=0.01)
        >>> p2 = rebound.Particle(simulation=sim, m=0.0, x=1., vy=1.)
        >>> p3 = rebound.Particle(simulation=sim, m=0.001, a=1.5, h=0.1, k=0.2, l=0.1)
        >>> p4 = rebound.Particle(simulation=sim, m=0.001, a=1.5, omega="uniform")  # omega will be a random number between 0 and 2pi

        """        

        # Unpickling
        binarydata = None
        if isinstance(simulation, bytes):
            # simulation is actually a bytes array with the particle data
            binarydata = simulation
        if isinstance(particle, bytes):
            # simulation is actually a bytes array with the particle data
            binarydata = particle
        if binarydata is not None:
            if len(binarydata) != sizeof(self):
                raise ValueError("Binary particle data does not have the right size.")
            buft = c_char * len(binarydata)
            buf = buft.from_buffer_copy(binarydata)
            memmove(byref(self), byref(buf), sizeof(self))
            self.c = 0
            self.sim = 0
            self.ap = 0
            return

        # Random initialization of particle angles
        clibrebound.reb_random_uniform.restype = c_double
        if simulation is not None:
            # Will use random seed stored in simulation.
            simp = byref(simulation)
        else:
            # Will use random seed based on time.
            simp = 0
        if Omega == "uniform":
            Omega = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if omega == "uniform":
            omega = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if pomega == "uniform":
            pomega = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if f == "uniform":
            f = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if M == "uniform":
            M = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if E == "uniform":
            E = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if l == "uniform":
            l = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if theta == "uniform":
            theta = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))
        if inc == "uniform":
            inc = clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(math.pi*2.0))

        self.hash = hash # set via the property, which checks for type
            
        if isinstance(primary, (str,int)):
           primary = simulation.particles[primary]


        if variation:
            if primary is None:
                primary = simulation.particles[0]
            # Find particle to differenciate
            lc = locals().copy()
            del lc["self"]
            del lc["variation"]
            del lc["variation2"]
            if particle is None:
                particle = Particle(**lc)
            # First or second order?
            if variation and variation2:
                variation_order = 2
            else:
                variation_order = 1
            # Shortcuts for variable names
            if variation == "l":
                variation = "lambda"
            if variation2 == "l":
                variation2 = "lambda"
            if variation == "i":
                variation = "inc"
            if variation2 == "i":
                variation2 = "inc"

            variationtypes = ["m","a","e","inc","omega","Omega","f","k","h","lambda","ix","iy"]
            if variation_order==1:
                if variation in variationtypes:
                    method = getattr(clibrebound, 'reb_particle_derivative_'+variation)
                    method.restype = Particle
                    p = method(c_double(simulation.G), primary, particle)
                else:
                    raise ValueError("Variational particles can only be initializes using the derivatives with respect to one of the following: %s."%", ".join(variationtypes))
            elif variation_order==2:
                if variation in variationtypes and variation2 in variationtypes:
                    # Swap variations if needed
                    vi1 = variationtypes.index(variation)
                    vi2 = variationtypes.index(variation2)
                    if vi2 < vi1:
                        variation, variation2 = variation2, variation

                    method = getattr(clibrebound, 'reb_particle_derivative_'+variation+'_'+variation2)
                    method.restype = Particle
                    p = method(c_double(simulation.G), primary, particle)
                else:
                    raise ValueError("Variational particles can only be initializes using the derivatives with respect to one of the following: %s."%", ".join(variationtypes))
            else:
                raise ValueError("Variational equations beyond second order are not implemented.")
            self.m = p.m
            self.x = p.x
            self.y = p.y
            self.z = p.z
            self.vx = p.vx
            self.vy = p.vy
            self.vz = p.vz
            return 

        if particle is not None:
            memmove(byref(self), byref(particle), sizeof(self))
            return
        cart = [x,y,z,vx,vy,vz]
        orbi = [primary,a,P,e,inc,Omega,omega,pomega,f,M,E,l,theta,T]
        if pal_h is not None:
            if h is not None:
                raise ValueError("Both h and pal_h were given.")
            h = pal_h
        if pal_k is not None:
            if k is not None:
                raise ValueError("Both k and pal_k were given.")
            k = pal_k
        if pal_ix is not None:
            if ix is not None:
                raise ValueError("Both ix and pal_ix were given.")
            ix = pal_ix
        if pal_iy is not None:
            if iy is not None:
                raise ValueError("Both iy and pal_iy were given.")
            iy = pal_iy
        pal  = [h,k,ix,iy]
       
        self.ax = 0.
        self.ay = 0.
        self.az = 0.
        if m is None:
            self.m = 0.
        else:
            self.m = m 
        if r is None:
            self.r = 0.
        else:
            self.r = r
        self.last_collision = 0.
        self.c = None
        self.ap = None
        
        if notNone([e,inc,omega,pomega,Omega,M,f,E,theta,T]) and notNone(pal):
            raise ValueError("You cannot mix Pal coordinates (h,k,ix,iy) with the following orbital elements: e,inc,Omega,omega,pomega,f,M,E,theta,T. If a longitude/anomaly is needed in Pal coordinates, use l.")
        if notNone(cart) and notNone(orbi):
                raise ValueError("You cannot pass cartesian coordinates and orbital elements (and/or primary) at the same time.")
        if notNone(orbi):
            if simulation is None:
                raise ValueError("Need to specify simulation when initializing particle with orbital elements.")
            if primary is None:
                clibrebound.reb_simulation_com.restype = Particle
                primary = clibrebound.reb_simulation_com(byref(simulation)) # this corresponds to adding in Jacobi coordinates
            if jacobi_masses is True:
                interior_mass = 0
                for p in simulation.particles:
                    interior_mass += p.m
                # orbit conversion uses mu=G*(p.m+primary.m) so set prim.m=Mjac-m so mu=G*Mjac
                primary.m = simulation.particles[0].m*(self.m + interior_mass)/interior_mass - self.m
            if a is None and P is None:
                raise ValueError("You need to pass either a semimajor axis or orbital period to initialize the particle using orbital elements.")
            if a is not None and P is not None:
                raise ValueError("You can pass either the semimajor axis or orbital period, but not both.")
            if a is None:
                a = (P**2*simulation.G*(primary.m + self.m)/(4.*math.pi**2))**(1./3.)
            if notNone(pal):
                # Pal orbital parameters
                if h is None:
                    h = 0.
                if k is None:
                    k = 0.
                if l is None:
                    l = 0.
                if ix is None:
                    ix = 0.
                if iy is None:
                    iy = 0.
                if((ix*ix + iy*iy) > 4.0):
                    raise ValueError("Passed (ix, iy) coordinates are not valid, squared sum exceeds 4.")
                clibrebound.reb_particle_from_pal.restype = Particle
                p = clibrebound.reb_particle_from_pal(c_double(simulation.G), primary, c_double(self.m), c_double(a), c_double(l), c_double(k), c_double(h), c_double(ix), c_double(iy))
            else:
                # Normal orbital parameters
                if e is None:
                    e = 0.
                if inc is None:
                    inc = 0.
                if Omega is None:               # we require that Omega be passed if you want to specify longitude of node
                    Omega = 0.

                pericenters = [omega, pomega]   # Need omega for C function. Can specify it either directly or through pomega indirectly. 
                numNones = pericenters.count(None)

                if numNones == 0:
                    raise ValueError("Can't pass both omega and pomega")
                if numNones == 2:                   # Neither passed.  Default to 0.
                    omega = 0.
                if numNones == 1:
                    if pomega is not None:          # Only have to find omega is pomega was passed
                        if math.cos(inc) > 0:       # inc is in range [-pi/2, pi/2] (prograde), so pomega = Omega + omega
                            omega = pomega - Omega
                        else:
                            omega = Omega - pomega  # for retrograde orbits, pomega = Omega - omega

                longitudes = [f,M,E,l,theta,T]      # can specify longitude through any of these.  Need f for C function.
                numNones = longitudes.count(None)

                if numNones < 5:
                    raise ValueError("Can only pass one longitude/anomaly in the set [f, M, E, l, theta, T]")
                if numNones == 6:                           # none of them passed.  Default to 0.
                    f = 0.
                if numNones == 5:                           # Only one was passed. Calculate f.
                    if f is not None:
                        pass                            # Nothing to be done
                    elif theta is not None:               # theta is next easiest
                        if math.cos(inc) > 0:           # for prograde orbits, theta = Omega + omega + f
                            f = theta - Omega - omega
                        else:
                            f = Omega - omega - theta   # for retrograde, theta = Omega - omega - f
                    elif l is not None:                 
                        if math.cos(inc) > 0:           # for prograde orbits, l = Omega + omega + M
                            M = l - Omega - omega
                        else:
                            M = Omega - omega - l       # for retrograde, l = Omega - omega - M
                        f = M_to_f(e, M)
                    elif T is not None:                 # works for both elliptical and hyperbolic orbits
                                                        # TODO: has accuracy problems for M=n*(t-T) << 1
                        n = (simulation.G*(primary.m+self.m)/abs(a**3))**0.5
                        M = n*(simulation.t - T)
                        f = M_to_f(e, M)
                    elif M is not None:
                        f = M_to_f(e, M)
                    elif E is not None:
                        f = E_to_f(e, E)

                err = c_int()
                clibrebound.reb_particle_from_orbit_err.restype = Particle
                p = clibrebound.reb_particle_from_orbit_err(c_double(simulation.G), primary, c_double(self.m), c_double(a), c_double(e), c_double(inc), c_double(Omega), c_double(omega), c_double(f), byref(err))
                if err.value == 1:
                    raise ValueError("Can't set e exactly to 1.")
                if err.value == 2:
                    raise ValueError("Eccentricity must be greater than or equal to zero.")
                if err.value == 3:
                    raise ValueError("Bound orbit (a > 0) must have e < 1.")
                if err.value == 4:
                    raise ValueError("Unbound orbit (a < 0) must have e > 1.")
                if err.value == 5:
                    raise ValueError("Unbound orbit can't have f beyond the range allowed by the asymptotes set by the hyperbola.")
                if err.value == 6:
                    raise ValueError("Primary has no mass.")
            self.x = p.x
            self.y = p.y
            self.z = p.z
            self.vx = p.vx
            self.vy = p.vy
            self.vz = p.vz
        else:
            if x is None:
                x = 0.
            if y is None:
                y = 0.
            if z is None:
                z = 0.
            if vx is None:
                vx = 0.
            if vy is None:
                vy = 0.
            if vz is None:
                vz = 0.
            self.x = x
            self.y = y
            self.z = z
            self.vx = vx
            self.vy = vy
            self.vz = vz
        
    def copy(self):
        """
        Returns a deep copy of the particle. The particle is not added to any simulation by default.
        """
        np = Particle()
        memmove(byref(np), byref(self), sizeof(self))
        return np

# Pickling method
    def __reduce__(self):
        return (Particle, (string_at(byref(self), size=sizeof(self)),))

    def orbit(self, primary=None, G=None):
        """ 
        Returns a rebound.Orbit object with the keplerian orbital elements
        corresponding to the particle around the passed primary
        (rebound.Particle) If no primary is passed, defaults to Jacobi coordinates
        (with mu = G*Minc, where Minc is the total mass from index 0 to the particle's index, inclusive). 
        
        Examples
        --------
        
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(x=1.,vy=1.)
        >>> orbit = sim.particles[1].orbit(sim.particles[0]) # Heliocentric coordinates
        >>> print(orbit.e) # gives the eccentricity

        Parameters
        ----------
        primary : rebound.Particle
            Central body (Optional. Default uses Jacobi coordinates)
        G : float
            Gravitational constant (Optional. Default takes G from simulation in which particle is in)
        
        Returns
        -------
        A rebound.Orbit object 
        """
        if not self._sim:
            # Particle not in a simulation
            if primary is None:
                raise ValueError("Particle does not belong to any simulation and no primary given. Cannot calculate orbit.")
            if G is None:
                raise ValueError("Particle does not belong to any simulation and G not given. Cannot calculate orbit.")
            else:
                G = c_double(G)
        else:
            # First check whether this is particles[0]
            clibrebound.reb_simulation_particle_index.restype = c_int
            index = clibrebound.reb_simulation_particle_index(byref(self)) # first check this isn't particles[0]
            if index == 0 and primary is None:
                raise ValueError("Orbital elements for particle[0] not implemented unless primary is provided")

            if primary is None:    # Use default, i.e., Jacobi coordinates
                clibrebound.reb_simulation_jacobi_com.restype = Particle   # now return jacobi center of mass
                primary = clibrebound.reb_simulation_jacobi_com(byref(self))
            G = c_double(self._sim.contents.G)
        
        err = c_int()
        clibrebound.reb_orbit_from_particle_err.restype = Orbit
        o = clibrebound.reb_orbit_from_particle_err(G, self, primary, byref(err))

        if err.value == 1:
            raise ValueError("Primary has no mass.")
        if err.value == 2:
            raise ValueError("Particle and primary positions are the same.")

        return o
    
    def sample_orbit(self, Npts=100, primary=None, samplingAngle=None, duplicateEndpoint=None):
        """
        Returns a nested list of xyz positions along the osculating orbit of the particle. 
        If primary is not passed, returns xyz positions along the Jacobi osculating orbit
        (with mu = G*Minc, where Minc is the total mass from index 0 to the particle's index, inclusive). 

        Parameters
        ----------
        Npts    : int, optional  
            Number of points along the orbit to return  (default: 100)
        primary : rebound.Particle, optional
            Primary to use for the osculating orbit (default: Jacobi center of mass)
        samplingAngle: str, optional
            This determines which angle is sampled linearly. Can be:
              - "M" (mean anomaly)
              - "E" (eccentric anomaly) 
              - "f" (true anomaly)
        duplicateEndpoint: bool, optional
            If true, then the first and last point will be identical for closed orbits. This is useful for some plotting tools. The default is true for eccentric orbits. The argument has no effect for hyperbolic orbits (because the endpoints are not identical).
        """
        if primary is None:
            primary = self.jacobi_com
        o = self.orbit(primary=primary)

        phases_f = []
        if samplingAngle is not None:
            if any(c not in "EMf" for c in samplingAngle):
                raise ValueError("Unknown character in samplingAngle.")
        
        if o.a < 0.: # hyperbolic orbit
            #a = o.a
            if samplingAngle is None:
                samplingAngle = "Mf"
            Nptsangle = {}
            for angle in samplingAngle[1:]:
                Nptsangle[angle] = (Npts-1)//len(samplingAngle) # one point is reserved for actual position
            Nptsangle[samplingAngle[0]] = Npts-1-sum(Nptsangle.values())
            if "M" in samplingAngle:
                phi = math.acos(-1./o.e)*0.999
                Npts = Nptsangle["M"]
                dphi = 2*phi/(Npts-1)
                for i in range(Npts):
                    f = M_to_f(o.e, phi)
                    phases_f.append(f)
                    phi -= dphi
            if "E" in samplingAngle:
                phi = math.acos(-1./o.e)*0.999
                Npts = Nptsangle["E"]
                dphi = 2*phi/(Npts-1)
                for i in range(Npts):
                    f = E_to_f(o.e, phi)
                    phases_f.append(f)
                    phi -= dphi
            if "f" in samplingAngle:
                phi = math.acos(-1./o.e)*0.999
                Npts = Nptsangle["f"]
                dphi = 2*phi/(Npts-1)
                for i in range(Npts):
                    f = mod2pi(phi)
                    phases_f.append(f)
                    phi -= dphi
        else:       # circular orbit
            #a = primary.m/(primary.m+self.m)*o.a
            if samplingAngle is None:
                samplingAngle = "Ef"
            if duplicateEndpoint is None:
                duplicateEndpoint = True
            Nptsangle = {}
            for angle in samplingAngle[1:]:
                Nptsangle[angle] = (Npts-1)//len(samplingAngle) # one point is reserved for actual position
            Nptsangle[samplingAngle[0]] = Npts-1-sum(Nptsangle.values())
            if "M" in samplingAngle:
                Npts = Nptsangle["M"]
                dphi = 2.*math.pi/(Npts-1 if duplicateEndpoint else Npts)  # one point is reserved for the end point
                for i in range(Npts):
                    f = M_to_f(o.e, i*dphi)
                    phases_f.append(f)
            if "E" in samplingAngle:
                Npts = Nptsangle["E"]
                dphi = 2.*math.pi/(Npts-1 if duplicateEndpoint else Npts)  # one point is reserved for the end point
                for i in range(Npts):
                    f = E_to_f(o.e, i*dphi)
                    phases_f.append(f)
            if "f" in samplingAngle:
                Npts = Nptsangle["f"]
                dphi = 2.*math.pi/(Npts-1 if duplicateEndpoint else Npts)  # one point is reserved for the end point
                for i in range(Npts):
                    f = i*dphi
                    f = mod2pi(f)
                    phases_f.append(f)

        # add actual position
        f = mod2pi(o.f)
        phases_f.append(f)
        phases_f.sort()
      
        pts_pre = []
        pts_post = []
        
        #clibrebound.reb_particle_com_of_pair.restype = Particle
        #primary = clibrebound.reb_particle_com_of_pair(primary, self)

        for f in phases_f:
            newp = Particle(a=o.a, f=f, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=self.m, primary=primary, simulation=self._sim.contents)
            if f<=o.f:
                pts_pre.append(newp.xyz)
            else:
                pts_post.append(newp.xyz)
        
        return pts_post + pts_pre

    # Simple operators for particles.
    
    def __eq__(self, other):
        # This ignores the pointer values
        if not isinstance(other,Particle):
            return NotImplemented
        clibrebound.reb_particle_diff.restype = c_int
        ret = clibrebound.reb_particle_diff(self, other)
        return not ret
    
    def __pow__(self, other):
        if not isinstance(other, Particle):
            return NotImplemented 
        clibrebound.reb_particle_distance.restype = c_double
        return clibrebound.reb_particle_distance(byref(self), byref(other))
    def __add__(self, other):
        if not isinstance(other, Particle):
            return NotImplemented 
        c = self.copy()
        return c.__iadd__(other)
    
    def __iadd__(self, other):
        if not isinstance(other, Particle):
            return NotImplemented 
        clibrebound.reb_particle_iadd(byref(self), byref(other))
        return self
    
    def __sub__(self, other):
        if not isinstance(other, Particle):
            return NotImplemented 
        c = self.copy()
        return c.__isub__(other)
    
    def __isub__(self, other):
        if not isinstance(other, Particle):
            return NotImplemented 
        clibrebound.reb_particle_isub(byref(self), byref(other))
        return self
    
    def __mul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented 
        c = self.copy()
        return c.__imul__(other)
    
    def __imul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented 
        clibrebound.reb_particle_imul(byref(self), c_double(other))
        return self
    
    def __rmul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented 
        c = self.copy()
        return c.__imul__(other)
    
    def __div__(self, other):
        return self.__truediv__(other)
    
    def __idiv__(self, other):
        return self.__itruediv__(other)
    
    def __truediv__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented 
        c = self.copy()
        if other==0.:
            raise ZeroDivisionError
        return c.__imul__(1./other)
    
    def __itruediv__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented 
        if other==0.:
            raise ZeroDivisionError
        return self.__imul__(1./other)
    
    def rotate(self, q):
        if not isinstance(q, Rotation):
            raise NotImplementedError
        clibrebound.reb_particle_irotate(byref(self), q)

    @property
    def index(self):
        clibrebound.reb_simulation_particle_index.restype = c_int
        return clibrebound.reb_simulation_particle_index(byref(self)) 
    
    @property
    def xyz(self):
        """
        Get or set the xyz position coordinates of the particle.
        """
        return [self.x, self.y, self.z]
    @xyz.setter
    def xyz(self, value):
        if len(value)!=3:
            raise AttributeError("Can only set xyz positions to array with length 3")
        self.x = float(value[0])
        self.y = float(value[1])
        self.z = float(value[2])
    
    @property
    def vxyz(self):
        """
        Get or set the xyz velocity coordinates of the particle.
        """
        return [self.vx, self.vy, self.vz]
    @vxyz.setter
    def vxyz(self, value):
        if len(value)!=3:
            raise AttributeError("Can only set xyz velocities to array with length 3")
        self.vx = float(value[0])
        self.vy = float(value[1])
        self.vz = float(value[2])
        
    @property
    def d(self):
        return self.orbit().d
    @property
    def v(self):
        return self.orbit().v 
    @property
    def h(self):
        return self.orbit().h
    @property
    def hvec(self):
        h = self.orbit().hvec
        return [h.x, h.y, h.z]
    @property
    def P(self):
        return self.orbit().P
    @P.setter
    def P(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, P=value, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def n(self):
        return self.orbit().n 
    @property
    def a(self):
        return self.orbit().a 
    @a.setter
    def a(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=value, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def rhill(self):
        return self.orbit().rhill
    @property
    def e(self):
        return self.orbit().e 
    @e.setter
    def e(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=value, inc=o.inc, omega=o.omega, Omega=o.Omega, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def evec(self):
        e = self.orbit().evec
        return [e.x, e.y, e.z]
    @property
    def inc(self):
        return self.orbit().inc 
    @inc.setter
    def inc(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=value, omega=o.omega, Omega=o.Omega, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def Omega(self):
        return self.orbit().Omega 
    @Omega.setter
    def Omega(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=o.omega, Omega=value, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def omega(self):
        return self.orbit().omega 
    @omega.setter
    def omega(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=value, Omega=o.Omega, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def pomega(self):
        return self.orbit().pomega 
    @pomega.setter
    def pomega(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, pomega=value, Omega=o.Omega, f=o.f) 
        self._copy_coordinates(newP)
    @property
    def f(self):
        return self.orbit().f 
    @f.setter
    def f(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, f=value) 
        self._copy_coordinates(newP)
    @property
    def M(self):
        return self.orbit().M 
    @M.setter
    def M(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, M=value) 
        self._copy_coordinates(newP)
    @property
    def l(self):
        return self.orbit().l 
    @l.setter
    def l(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, l=value) 
        self._copy_coordinates(newP)
    @property
    def theta(self):
        return self.orbit().theta 
    @theta.setter
    def theta(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, theta=value) 
        self._copy_coordinates(newP)
    @property
    def T(self):
        return self.orbit().T
    @T.setter
    def T(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, e=o.e, inc=o.inc, omega=o.omega, Omega=o.Omega, T=value) 
        self._copy_coordinates(newP)
    # Pal coordinates
    @property
    def pal_h(self):
        return self.orbit().pal_h
    @pal_h.setter
    def pal_h(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, l=o.l, pal_h=value, pal_k=o.pal_k, pal_ix=o.pal_ix, pal_iy=o.pal_iy)
        self._copy_coordinates(newP)
    @property
    def pal_k(self):
        return self.orbit().pal_k
    @pal_k.setter
    def pal_k(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, l=o.l, pal_h=o.pal_h, pal_k=value, pal_ix=o.pal_ix, pal_iy=o.pal_iy)
        self._copy_coordinates(newP)
    @property
    def pal_ix(self):
        return self.orbit().pal_ix
    @pal_ix.setter
    def pal_ix(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, l=o.l, pal_h=o.pal_h, pal_k=o.pal_k, pal_ix=value, pal_iy=o.pal_iy)
        self._copy_coordinates(newP)
    @property
    def pal_iy(self):
        return self.orbit().pal_iy
    @pal_iy.setter
    def pal_iy(self,value):
        o = self.orbit()
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        primary = clibrebound.reb_simulation_jacobi_com(byref(self))
        if self._sim is None:
            raise RuntimeError("Cannot modify particle which is not a member of a simulation.")
        newP = Particle(simulation=self._sim.contents, primary=primary, m=self.m, a=o.a, l=o.l, pal_h=o.pal_h, pal_k=o.pal_k, pal_ix=o.pal_ix, pal_iy=value)
        self._copy_coordinates(newP)
    # Other properties
    @property
    def jacobi_com(self):
        clibrebound.reb_simulation_jacobi_com.restype = Particle
        return clibrebound.reb_simulation_jacobi_com(byref(self))
    @property
    def hash(self):
        """
        Get or set the particle's hash.  If set to a string, the corresponding integer hash is calculated.
        """
        return c_uint32(self._hash)
    @hash.setter
    def hash(self, value):
        PY3 = sys.version_info[0] == 3
        hash_types = c_uint32, c_uint, c_uint64
        if PY3:
            string_types = str,
            int_types = int,
        else:
            string_types = basestring,
            int_types = int, long,
        if isinstance(value, hash_types):
            self._hash = value.value
        elif isinstance(value, string_types):
            self._hash = hash(value).value
        elif isinstance(value, int_types):
            self._hash = value
        else:
            raise AttributeError("Hash must be set to an integer, a ctypes.c_uint32 or a string. See UniquelyIdentifyingParticlesWithHashes.ipynb ipython_example.")

    def _copy_coordinates(self, p):
        """
        Copy coordinates (and only coordinates) from particle p to self
        """
        self.xyz = p.xyz
        self.vxyz = p.vxyz

from .simulation import Simulation
from . import clibrebound
from .tools import E_to_f, M_to_f, mod2pi
from .orbit import Orbit
from .rotation import Rotation
from .hash import hash

if sizeof(c_void_p)==4:
    # Add padding for 32 bit
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
                    ("last_collision", c_double),
                    ("c", c_void_p),
                    ("_pad1", c_char*4),
                    ("_hash", c_uint32),
                    ("_pad2", c_char*4),
                    ("ap", c_void_p),
                    ("_pad3", c_char*4),
                    ("_sim", POINTER(Simulation)),
                    ("_pad4", c_char*4),
                         ]
else:
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
                    ("last_collision", c_double),
                    ("c", c_void_p),
                    ("_hash", c_uint32),
                    ("ap", c_void_p),
                    ("_sim", POINTER(Simulation)),
                         ]
