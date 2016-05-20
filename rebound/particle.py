from ctypes import Structure, c_double, c_int, byref, memmove, sizeof
from . import clibrebound
import math
import ctypes.util
import rebound

__all__ = ["Particle"]

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
    lastcollision : float       
        Last time the particle had a physical collision (if checking for collisions)
    c           : c_void_p (C void pointer) 
        Pointer to the cell the particle is currently in (if using tree code)
    id          : int         
        Particle ID (arbitrary, specified by the user)
    ap          : c_void_p (C void pointer)
        Pointer to additional parameters one might want to add to particles
    _sim        : POINTER(rebound.Simulation)
        Internal pointer to the parent simulation (used in C version of REBOUND)
    """
    def __str__(self):
        """ 
        Returns a string with the position and velocity of the particle.
        """
        return "<rebound.Particle object, id=%s m=%s x=%s y=%s z=%s vx=%s vy=%s vz=%s>"%(self.id,self.m,self.x,self.y,self.z,self.vx,self.vy,self.vz)
    
    def __init__(self, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, P=None, e=None, inc=None, Omega=None, omega=None, pomega=None, f=None, M=None, l=None, theta=None, T=None, r=None, id=None, date=None, simulation=None, variation=None, variation2=None, h=None, k=None, ix=None, iy=None):
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
        
        For classical orbital paramerers, one can specify the longitude 
        of the ascending node by passing Omega, to specify the pericenter 
        one can pass either omega or pomega (not both), and for the 
        longitude/anomaly one can pass one of f, M, l or theta.  
        See ipython_examples/OrbitalElements.ipynb for examples.  
        See also Murray & Dermott Solar System Dynamics for formal 
        definitions of angles in orbital mechanics.

        All angles should be specified in radians.
        

        Parameters
        ----------
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
        l           : float       
            Mean longitude              (Default: 0)
        theta       : float       
            True longitude              (Default: 0)
        T           : float 
            Time of pericenter passage  
        h           : float       
            h variable, see Pal (2009) for a definition  (Default: 0)
        k           : float       
            k variable, see Pal (2009) for a definition  (Default: 0)
        ix          : float       
            ix variable, see Pal (2009) for a definition  (Default: 0)
        iy           : float       
            iy variable, see Pal (2009) for a definition  (Default: 0)
        r           : float       
            Particle radius (only used for collisional simulations)
        id          : int               (Default: 0)
            Particle ID (arbitrary, specified by the user)
        date        : string      
            For consistency with adding particles through horizons.  Not used here.
        simulation  : Simulation  
            Simulation instance associated with this particle (Required)
        variation   : string            (Default: None)
            Set this string to the name of an orbital parameter to initialize the particle as a variational particle.
            Can be one of the following: m, a, e, inc, omega, Omega, f, k, h, lambda, ix, iy.
        variation2  : string            (Default: None)
            Set this string to the name of a second orbital parameter to initialize the particle as a second order variational particle. Only used for second order variational equations. 
            Can be one of the following: m, a, e, inc, omega, Omega, f, k, h, lambda, ix, iy.
        
        Returns
        -------
        A rebound.Particle object 
        
        Examples
        --------

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> p1 = rebound.Particle(simulation=sim, m=0.001, a=0.5, e=0.01)
        >>> p2 = rebound.Particle(simulation=sim, m=0.0, x=1., vy=1.)
        >>> p3 = rebound.Particle(simulation=sim, m=0.001, a=1.5, h=0.1, k=0.2, l=0.1)

        """        
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
                    method = getattr(clibrebound, 'reb_derivatives_'+variation)
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

                    method = getattr(clibrebound, 'reb_derivatives_'+variation+'_'+variation2)
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
        orbi = [primary,a,P,e,inc,Omega,omega,pomega,f,M,l,theta,T]
        pal  = [h,k,ix,iy]
       
        self.ax = 0.
        self.ay = 0.
        self.az = 0.
        if m is None:
            self.m = 0.
        else:
            self.m = m
        if id is None:
            self.id = 0
        else:
            self.id = id
        if r is None:
            self.r = 0.
        else:
            self.r = r
        self.lastcollision = 0.
        self.c = None
        self.ap = None
        
        if notNone(cart) and notNone(orbi):
                raise ValueError("You cannot pass cartesian coordinates and orbital elements (and/or primary) at the same time.")
        if notNone(orbi):
            if simulation is None:
                raise ValueError("Need to specify simulation when initializing particle with orbital elements.")
            if primary is None:
                clibrebound.reb_get_com.restype = Particle
                primary = clibrebound.reb_get_com(byref(simulation)) # this corresponds to adding in Jacobi coordinates
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
                clibrebound.reb_tools_pal_to_particle.restype = Particle
                p = clibrebound.reb_tools_pal_to_particle(c_double(simulation.G), primary, c_double(self.m), c_double(a), c_double(l), c_double(k), c_double(h), c_double(ix), c_double(iy))
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

                longitudes = [f,M,l,theta,T]      # can specify longitude through any of these four.  Need f for C function.
                numNones = longitudes.count(None)

                if numNones < 4:
                    raise ValueError("Can only pass one longitude/anomaly in the set [f, M, l, theta, T]")
                if numNones == 5:                           # none of them passed.  Default to 0.
                    f = 0.
                if numNones == 4:                           # Only one was passed.
                    if f is None:                           # Only have to work if f wasn't passed.
                        if theta is not None:               # theta is next easiest
                            if math.cos(inc) > 0:           # for prograde orbits, theta = Omega + omega + f
                                f = theta - Omega - omega
                            else:
                                f = Omega - omega - theta   # for retrograde, theta = Omega - omega - f
                        else:                               # Either M, l, or T was passed.  Will need to find M first (if not passed) to find f
                            if l is not None:
                                if math.cos(inc) > 0:       # for prograde orbits, l = Omega + omega + M
                                    M = l - Omega - omega
                                else:
                                    M = Omega - omega - l   # for retrograde, l = Omega - omega - M
                            else:
                                if T is not None:           # works for both elliptical and hyperbolic orbits
                                    n = (simulation.G*(primary.m+self.m)/abs(a**3))**0.5
                                    M = n*(simulation.t - T)
                            clibrebound.reb_tools_M_to_f.restype = c_double
                            f = clibrebound.reb_tools_M_to_f(c_double(e), c_double(M))

                err = c_int()
                clibrebound.reb_tools_orbit_to_particle_err.restype = Particle
                p = clibrebound.reb_tools_orbit_to_particle_err(c_double(simulation.G), primary, c_double(self.m), c_double(a), c_double(e), c_double(inc), c_double(Omega), c_double(omega), c_double(f), byref(err))
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
            self.m = p.m
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

    def calculate_orbit(self, primary=None, G=None):
        """ 
        Returns a rebound.Orbit object with the keplerian orbital elements
        corresponding to the particle around the passed primary
        (rebound.Particle) If no primary is passed, defaults to Jacobi coordinates. 
        
        Examples
        --------
        
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(x=1.,vy=1.)
        >>> orbit = sim.particles[1].calculate_orbit(sim.particles[0])
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
            clibrebound.reb_get_particle_index.restype = c_int
            index = clibrebound.reb_get_particle_index(byref(self)) # first check this isn't particles[0]
            if index == 0:
                raise ValueError("Orbital elements for particle[0] not implemented.")

            if primary is None:    # Use default, i.e., Jacobi coordinates
                clibrebound.reb_get_jacobi_com.restype = Particle   # now return jacobi center of mass
                primary = clibrebound.reb_get_jacobi_com(byref(self))
            G = c_double(self._sim.contents.G)
        
        err = c_int()
        clibrebound.reb_tools_particle_to_orbit_err.restype = rebound.Orbit
        o = clibrebound.reb_tools_particle_to_orbit_err(G, self, primary, byref(err))

        if err.value == 1:
            raise ValueError("Primary has no mass.")
        if err.value == 2:
            raise ValueError("Particle and primary positions are the same.")

        return o
    
    # Simple operators for particles.

    def __sub__(self, other):
        if isinstance(other, Particle):
            clibrebound.reb_particle_minus.restype = Particle
            return clibrebound.reb_particle_minus(self, other)
        return NotImplemented 
    
    def __add__(self, other):
        if isinstance(other, Particle):
            clibrebound.reb_particle_plus.restype = Particle
            return clibrebound.reb_particle_plus(self, other)
        return NotImplemented 
    
    def __mul__(self, other):
        if isinstance(other, float):
            clibrebound.reb_particle_multiply.restype = Particle
            return clibrebound.reb_particle_multiply(self, c_double(other))
        return NotImplemented 
    
    def __rmul__(self, other):
        if isinstance(other, float):
            clibrebound.reb_particle_multiply.restype = Particle
            return clibrebound.reb_particle_multiply(self, c_double(other))
        return NotImplemented 

    def __truediv__(self, other):
        return self.__div__(other)

    def __div__(self, other):
        if isinstance(other, float):
            clibrebound.reb_particle_divide.restype = Particle
            return clibrebound.reb_particle_divide(self, c_double(other))
        return NotImplemented 

    @property
    def index(self):
        clibrebound.reb_get_particle_index.restype = c_int
        return clibrebound.reb_get_particle_index(byref(self)) 
    
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
        return self.calculate_orbit().d
    @property
    def v(self):
        return self.calculate_orbit().v 
    @property
    def h(self):
        return self.calculate_orbit().h
    @property
    def P(self):
        return self.calculate_orbit().P
    @property
    def n(self):
        return self.calculate_orbit().n 
    @property
    def a(self):
        return self.calculate_orbit().a 
    @property
    def e(self):
        return self.calculate_orbit().e 
    @property
    def inc(self):
        return self.calculate_orbit().inc 
    @property
    def Omega(self):
        return self.calculate_orbit().Omega 
    @property
    def omega(self):
        return self.calculate_orbit().omega 
    @property
    def pomega(self):
        return self.calculate_orbit().pomega 
    @property
    def f(self):
        return self.calculate_orbit().f 
    @property
    def M(self):
        return self.calculate_orbit().M 
    @property
    def l(self):
        return self.calculate_orbit().l 
    @property
    def theta(self):
        return self.calculate_orbit().theta 
    @property
    def T(self):
        return self.calculate_orbit().T
    @property
    def orbit(self):
        return self.calculate_orbit()

    @property
    def jacobi_com(self):
        clibrebound.reb_get_jacobi_com.restype = Particle
        return clibrebound.reb_get_jacobi_com(byref(self))
