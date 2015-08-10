from ctypes import *
from . import clibrebound
import math
import ctypes.util
import rebound

__all__ = ["Orbit", "Particle"]
TINY=1.e-308
MIN_REL_ERROR = 1.e-12

# Helper functions
TWOPI = 2.*math.pi
def mod2pi(f):
    """Returns the angle f (in radians) modulo 2 pi."""
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def eccentricAnomaly(e,M):
    """Returns the eccentric anomaly given the eccentricity and mean anomaly (in radians) of a Keplerian orbit.
    
    Parameters
    ----------
    e       : (float)       Eccentricity
    M       : (float)       Mean anomaly (in radians)
    """
    if e<1.:
        E = M if e<0.8 else math.pi
        
        F = E - e*math.sin(M) - M
        for i in range(100):
            E = E - F/(1.0-e*math.cos(E))
            F = E - e*math.sin(E) - M
            if math.fabs(F)<1e-16:
                break
        E = mod2pi(E)
        return E
    else:
        E = M 
        
        F = E - e*math.sinh(E) - M
        for i in range(100):
            E = E - F/(1.0-e*math.cosh(E))
            F = E - e*math.sinh(E) - M
            if math.fabs(F)<1e-16:
                break
        E = mod2pi(E)
        return E



def notNone(a):
    """Returns True if array a contains at least one element that is not None. Returns False otherwise."""
    return a.count(None) != len(a)


class Orbit():
    """
    A class containing orbital parameters for a particle.
    This is an abstraction of the reb_orbit data structure in C.

    When using the various REBOUND functions using Orbits, all angles are in radians. 

    Parameters
    ---------
    a       : (float)           semimajor axis
    r       : (float)           radial distance from reference 
    h       : (float)           specific angular momentum
    P       : (float)           orbital period
    l       : (float)           mean longitude = Omega + omega + M
    e       : (float)           eccentricity
    inc     : (float)           inclination
    Omega   : (float)           longitude of ascending node
    omega   : (float)           argument of pericenter
    f       : (float)           true anomaly
    """
    def __init__(self):
        self.a      =   None    # semimajor axis
        self.r      =   None    # radial distance from reference
        self.h      =   None    # angular momentum
        self.P      =   None    # orbital period
        self.l      =   None    # mean longitude = Omega + omega + M
        self.e      =   None    # eccentricity
        self.inc    =   None    # inclination
        self.Omega  =   None    # longitude of ascending node
        self.omega  =   None    # argument of perihelion
        self.f      =   None    # true anomaly

    def __str__(self):
        """
        Returns a string with the semi-major axis and eccentricity of the orbit.
        """
        return "<rebound.Orbit instance, a=%s e=%s>"%(str(self.a),str(self.e))


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
    
    def __init__(self, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, r=None, id=None, date=None, simulation=None):   
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
        primary     : (Particle)    Primary body for converting orbital elements to cartesian (Default: center of mass of the particles in the passed simulation) 
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
            raise ValueError("Cannot initialise particle from other particles.")
        cart = [x,y,z,vx,vy,vz]
        orbi = [primary,a,anom,e,omega,inc,Omega,MEAN]
        if m is None:   #default value for mass
            m = 0.
        if notNone(cart) and notNone(orbi):
                raise ValueError("You cannot pass cartesian coordinates and orbital elements (and/or primary) at the same time.")
        if notNone(orbi):
            if simulation is None:
                raise ValueError("Need to specify simulation when initializing particle with orbital elements.")
            if primary is None:
                clibrebound.reb_get_com.restype = Particle
                primary = clibrebound.reb_get_com(simulation.simulation)
            if a is None:
                raise ValueError("You need to pass a semi major axis to initialize the particle using orbital elements.")
            if anom is None:
                anom = 0.
            if e is None:
                e = 0.
            if omega is None:
                omega = 0.
            if inc is None:
                inc = 0.
            if Omega is None:
                Omega = 0.
            if MEAN is None:
                MEAN = False
            self.set_orbit(simulation, m=m,primary=primary,a=a,anom=anom,e=e,omega=omega,inc=inc,Omega=Omega,MEAN=MEAN)
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
            self.m = m
            self.x = x
            self.y = y
            self.z = z
            self.vx = vx
            self.vy = vy
            self.vz = vz

        self.id = -1 if id is None else id
        self.r  =  0. if r is None else r

    def set_orbit(self,
                    simulation,
                    m,          # mass
                    primary,    # central body (rebound.Particle object)
                    a,          # semimajor axis
                    anom=0.,    # anomaly
                    e=0.,       # eccentricity
                    omega=0.,   # argument of pericenter
                    inc=0.,     # inclination
                    Omega=0.,   # longitude of ascending node
                    MEAN=False):    # mean anomaly
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
    
        self.m = m

        if not(0.<=inc<=math.pi): raise ValueError('inc must be in range [0,pi]')
        if e>1.:
            if math.fabs(anom)>math.acos(-1./e): raise ValueError('hyperbolic orbit with anomaly larger than angle of asymptotes')
    
        if MEAN is True: # need to calculate f
            E = eccentricAnomaly(e,anom)
            if e>1.:
                print("not working yet")
                exit(1)
                f = 2.*math.atan(math.sqrt(-(1.+ e)/(1. - e))*math.tanh(0.5*E))
            else:
                f = mod2pi(2.*math.atan(math.sqrt((1.+ e)/(1. - e))*math.tan(0.5*E)))
        else:
            f = anom
        
        cO = math.cos(Omega)
        sO = math.sin(Omega)
        co = math.cos(omega)
        so = math.sin(omega)
        cf = math.cos(f)
        sf = math.sin(f)
        ci = math.cos(inc)
        si = math.sin(inc)
        
        r = a*(1.-e**2)/(1.+e*cf)
        
        # Murray & Dermott Eq. 2.122
        self.x  = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
        self.y  = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
        self.z  = primary.z + r*(so*cf+co*sf)*si
        
        n = math.sqrt(simulation.simulation.contents.G*(primary.m+m)/(a**3))
        if e>1.:
            v0 = n*a/math.sqrt(-(1.-e**2))
        else:
            v0 = n*a/math.sqrt(1.-e**2)
        
        # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
        self.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
        self.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
        self.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)


    def calculate_orbit(self, simulation, primary=None, index=None, verbose=False):
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
        index   : (int)              Index of particle in particles array (will return Jacobi elements) (Default: None)    
        verbose : (boolean)          If set to True, will print out error msgs (Default: False)
        
        Returns
        -------
        A rebound.Orbit object 
        """
        if primary is None:
            if index is None:
                primary = simulation.particles[0]
            else:
                primary = simulation.calculate_com(index)

        o = Orbit()
        if primary.m <= TINY:
            if verbose is True:
                print("Star has no mass.")
            return o                            # all values set to None
        
        dx = self.x - primary.x
        dy = self.y - primary.y
        dz = self.z - primary.z
        o.r = math.sqrt ( dx*dx + dy*dy + dz*dz )
        if o.r <= TINY:
            if verbose is True:
                print('Particle and primary positions are the same.')
            return o
        
        dvx = self.vx - primary.vx
        dvy = self.vy - primary.vy
        dvz = self.vz - primary.vz
        v = math.sqrt ( dvx*dvx + dvy*dvy + dvz*dvz )
        
        mu = simulation.simulation.contents.G*(self.m+primary.m)
        o.a = -mu/( v*v - 2.*mu/o.r )               # semi major axis
        
        h0 = (dy*dvz - dz*dvy)                      # angular momentum vector
        h1 = (dz*dvx - dx*dvz)
        h2 = (dx*dvy - dy*dvx)
        o.h = math.sqrt ( h0*h0 + h1*h1 + h2*h2 )   # abs value of angular momentum
        if o.h/(o.r*v) <= MIN_REL_ERROR:
            if verbose is True:
                print('Particle orbit is radial.')
            return o
        
        vr = (dx*dvx + dy*dvy + dz*dvz)/o.r
        e0 = 1./mu*( (v*v-mu/o.r)*dx - o.r*vr*dvx )
        e1 = 1./mu*( (v*v-mu/o.r)*dy - o.r*vr*dvy )
        e2 = 1./mu*( (v*v-mu/o.r)*dz - o.r*vr*dvz )
        o.e = math.sqrt( e0*e0 + e1*e1 + e2*e2 )   # eccentricity
        
        o.P = math.copysign(2.*math.pi*math.sqrt( math.fabs(o.a*o.a*o.a/mu) ), o.a)  # period
        o.inc = math.acos( h2/o.h )               # inclination (wrt xy-plane)
        # if pi/2<i<pi it's retrograde
        n0 = -h1                                    # node vector
        n1 =  h0                                    # in xy plane => no z component
        n = math.sqrt( n0*n0 + n1*n1 )
        er = dx*e0 + dy*e1 + dz*e2
        if n/(o.r*v)<=MIN_REL_ERROR or o.inc<=MIN_REL_ERROR:# we are in the xy plane
            o.Omega=0.
            if o.e <= MIN_REL_ERROR:              # omega not defined for circular orbit
                o.omega = 0.
            else:
                if e1>=0.:
                    o.omega=math.acos(e0/o.e)
                else:
                    o.omega = 2.*math.pi-math.acos(e0/o.e)
        else:
            if o.e <= MIN_REL_ERROR:
                o.omega = 0.
            else:
                if e2>=0.:                        # omega=0 if perictr at asc node
                    o.omega=math.acos(( n0*e0 + n1*e1 )/(n*o.e))
                else:
                    o.omega=2.*math.pi-math.acos(( n0*e0 + n1*e1 )/(n*o.e))
            if n1>=0.:
                o.Omega = math.acos(n0/n)
            else:
                o.Omega=2.*math.pi-math.acos(n0/n)# Omega=longitude of asc node
        # taken in xy plane from x axis
        
        if o.e<=MIN_REL_ERROR:                           # circular orbit
            o.f=0.                                  # f has no meaning
            o.l=0.
        else:
            cosf = er/(o.e*o.r)
            cosea = (1.-o.r/o.a)/o.e
            
            if -1.<=cosf and cosf<=1.:                       # failsafe
                o.f = math.acos(cosf)
            else:
                o.f = math.pi/2.*(1.-cosf)
            
            if -1.<=cosea and cosea<=1.:
                ea  = math.acos(cosea)
            else:
                ea = math.pi/2.*(1.-cosea)
            
            if vr<0.:
                o.f=2.*math.pi-o.f
                ea =2.*math.pi-ea
            
            o.l = ea -o.e*math.sin(ea) + o.omega+ o.Omega  # mean longitude
        
        return o

