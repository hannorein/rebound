from ctypes import *
from . import clibrebound
import math
import ctypes.util
import rebound

__all__ = ["Particle"]

def notNone(a):
    """Returns True if array a contains at least one element that is not None. Returns False otherwise."""
    return a.count(None) != len(a)

class Particle(Structure):
    """
    The main REBOUND particle data structure. 
    This is an abstraction of the reb_particle structure in C
    
    Parameters
    ----------
    x, y, z     : (float)       Particle positions
    vx, vy, vz  : (float)       Particle velocities
    ax, ay, az  : (float)       Particle accelerations
    m           : (float)       Particle mass
    r           : (float)       Particle radius
    lastcollision:(float)       Last time the particle had a physical collision (if checking for collisions)
    c           : (float)       Pointer to the cell the particle is currently in (if using tree code)
    id          : (int)         Particle ID (arbitrary, specified by the user)
    """
    _fields_ = [("x", c_double),
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
                ("id", c_int)]
    def __str__(self):
        """ 
        Returns a string with the position and velocity of the particle.
        """
        return "<rebound.Particle object, id=%s m=%s x=%s y=%s z=%s vx=%s vy=%s vz=%s>"%(self.id,self.m,self.x,self.y,self.z,self.vx,self.vy,self.vz)
    
    def __init__(self, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, e=None, inc=None, Omega=None, omega=None, f=None, pomega=None, M=None, l=None, theta=None, r=None, id=None, date=None, simulation=None):
        """
        Initializes a Particle structure.
        Typically users will not create Particle structures directly.
        Rather, use the add member function of a Simulation instance, which will create a Particle and add it to the simulation.

        Accepts either cartesian positions and velocities, or orbital elements together with the reference Particle (primary), but not both.
        Requires a simulation instance for the Particle.
        For convenience, optional keywords that are not passed default to zero (mass, cartesian and orbital elements). 
        However, if passing orbital elements directly or implicitly (by passing a primary Particle), you must specify the semimajor axis.
        All angles should be specified in radians.

        Parameters
        ---------
        particle    : (Particle)    For consistency with other particle addition routines.  Cannot be passed when creating a particle in this way.
        m           : (float)       Mass        (Default: 0)
        x, y, z     : (float)       Positions   (Default: 0)
        vx, vy, vz  : (float)       Velocities  (Default: 0)
        primary     : (Particle)    Primary body for converting orbital elements to cartesian (Default: center of mass of the particles in the passed simulation, i.e., will yield Jacobi coordinates as you progressively pass particles) 
        a           : (float)       Semimajor axis (Required if passing orbital elements)
        anom        : (float)       Either the true or mean anomaly, determined by the value of the MEAN keyword (Default: 0)
        e           : (float)       Eccentricity                (Default: 0)
        omega       : (float)       Argument of pericenter      (Default: 0)
        inc         : (float)       Inclination                 (Default: 0)
        Omega       : (float)       Longitude of ascending node (Default: 0)
        MEAN        : (boolean)     Flag for whether anom refers to the mean anomaly (MEAN=True) or true anomaly (MEAN=False) (Default: False)
        r           : (float)       Particle radius (only used for collisional simulations)
        id          : (int)         Particle ID (arbitrary, specified by the user)
        date        : (string)      For consistency with adding particles through horizons.  Not used here.
        simulation  : (Simulation)  Simulation instance associated with this particle (Required)
        
        Usage
        -----
        sim = rebound.Simulation()
        primary = rebound.Particle(m=1., simulation=sim) # particle with unit mass at origin & v=0
        
        # test particle (m=0) with specified elements using mean anomaly
        p = rebound.Particle(m=0.,primary=primary,a=2.5, anom=math.pi/2,e=0.3,omega=math.pi/6,inc=math.pi/3,Omega=0.,MEAN=True, simulation=sim)
        
        # m=0.1 particle on circular orbit on positive x axis 
        p2 = rebound.Particle(m=0.1, x=1., vy=1., simulation=sim)
        """        

        if particle is not None:
            raise ValueError("Cannot initialize particle from other particles.")
        cart = [x,y,z,vx,vy,vz]
        orbi = [primary,a,e,inc,Omega,omega,pomega,f,M,l,theta]
        
        self.m = 0. if m is None else m
        self.id = -1 if id is None else id
        self.r  =  0. if r is None else r
        
        if notNone(cart) and notNone(orbi):
                raise ValueError("You cannot pass cartesian coordinates and orbital elements (and/or primary) at the same time.")
        if notNone(orbi):
            if simulation is None:
                raise ValueError("Need to specify simulation when initializing particle with orbital elements.")
            if primary is None:
                clibrebound.reb_get_com.restype = Particle
                primary = clibrebound.reb_get_com(byref(simulation))
            if a is None:
                raise ValueError("You need to pass a semi major axis to initialize the particle using orbital elements.")
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

            longitudes = [f,M,l,theta]      # can specify longitude through any of these four.  Need f for C function.
            numNones = longitudes.count(None)

            if numNones < 3:
                raise ValueError("Can only pass one longitude/anomaly in the set [f, M, l, theta]")
            if numNones == 4:                           # none of them passed.  Default to 0.
                f = 0.
            if numNones == 3:                           # Only one was passed.
                if f is None:                           # Only have to work if f wasn't passed.
                    if theta is not None:               # theta is next easiest
                        if math.cos(inc) > 0:           # for prograde orbits, theta = Omega + omega + f
                            f = theta - Omega - omega
                        else:
                            f = Omega - omega - theta   # for retrograde, theta = Omega - omega - f
                    else:                               # Either M or l was passed.  Will need to find M first to find f.
                        if l is not None:
                            if math.cos(inc) > 0:       # for prograde orbits, l = Omega + omega + M
                                M = l - Omega - omega
                            else:
                                M = Omega - omega - l   # for retrograde, l = Omega - omega - M
                        clibrebound.reb_tools_M_to_f.restype = c_double
                        f = clibrebound.reb_tools_M_to_f(c_double(e), c_double(M)).value

            err = c_int()
            clibrebound.reb_tools_orbit_to_particle_err.restype = Particle
            p = clibrebound.reb_tools_orbit_to_particle_err(c_double(simulation.G), primary, c_double(self.m), c_double(a), c_double(e), c_double(inc), c_double(Omega), c_double(omega), c_double(f), byref(err))

            if err.value == 1:
                raise ValueError("Can't initialize a radial orbit with orbital elements.")
            if err.value == 2:
                raise ValueError("Eccentricity must be greater than or equal to zero.")
            if err.value == 3:
                raise ValueError("Bound orbit (a > 0) must have e < 1.")
            if err.value == 4:
                raise ValueError("Unbound orbit (a < 0) must have e > 1.")
            if err.value == 5:
                raise ValueError("Unbound orbit can't have f beyond the range allowed by the asymptotes set by the hyperbola.")
            
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


    def set_orbit(self,
                    simulation,
                    m,          # mass
                    primary,    # central body (rebound.Particle object)
                    a,          # semimajor axis
                    e=0.,       # eccentricity
                    inc=0.,     # inclination
                    Omega=0.,   # longitude of ascending node
                    omega=0.,   # argument of pericenter
                    f=0.):      # true anomaly
        """ 
        Initialises a particle structure with the passed set of
        orbital elements. Mass (m), primary and 'a' are required (see Parameters
        below, and any orbital mechanics text, e.g., Murray & Dermott
        Solar System Dynamics for definitions). Other values default to zero.
        All angles should be passed in radians. Units are set by the
        gravitational constant G (default = 1.). If MEAN is set to True, anom is
        taken as the mean anomaly, rather than the true anomaly.
        
        Parameters
        ----------
        m       : (float)           Mass of the particle (Required)
        primary : (rebound.Particle)Particle structure for the central body (Required)
        a       : (float)           Semimajor axis (Required)
        anom    : (float)           Either the true or mean anomaly, determined by the value of the MEAN keyword (Default: 0)
        e       : (float)           Eccentricity (Default: 0)
        omega   : (float)           Argument of pericenter (Default: 0)
        inc     : (float)           Inclination (to xy plane) (Default: 0)
        Omega   : (float)           Longitude of the ascending node (Default: 0)
        MEAN    : (boolean)         If False, anom = true anomaly. If True, anom = mean anomaly (Default: False)
        
        Returns
        -------
        A rebound.Particle structure initialized with the given orbital parameters
        """
   
    def calculate_orbit(self, simulation, primary=None):
        """ 
        Returns a rebound.Orbit object with the keplerian orbital elements
        corresponding to the particle around the central body primary
        (rebound.Particle). By default will return orbital elements referenced
        to simulation.particles[0].  If an index is passed instead of primary,
        it will return Jacobi orbital elements
        (with primary as the center of mass of all interior particles).
        Edge cases will return values set to None. If
        verbose is set to True (default=False), error messages are printed
        when a breakout condition is met.
        
        Usage
        -----
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(x=1.,vy=1.)

        orbit = sim.particles[1].calculate_orbit(sim)
        print(orbit.e) # gives the eccentricity

        orbit = sim.particles[1].calculate_orbit(sim, verbose=True)

        Parameters
        ----------
        primary : (rebound.Particle) Central body (Default: COM of interior bodies)
        verbose : (boolean)          If set to True, will print out error msgs (Default: False)
        
        Returns
        -------
        A rebound.Orbit object 
        """
        if primary is None:
            primary = simulation.particles[0]
        
        err = c_int()
        clibrebound.reb_tools_particle_to_orbit_err.restype = rebound.Orbit
        o = clibrebound.reb_tools_particle_to_orbit_err(c_double(simulation.G), self, primary, byref(err))

        if err.value == 1:
            raise ValueError("Primary has no mass.")
        if err.value == 2:
            raise ValueError("Particle and primary positions are the same.")

        return o
