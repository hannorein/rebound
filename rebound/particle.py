from ctypes import *
import math
import ctypes.util
import rebound

__all__ = ["Orbit", "Particle"]
TINY=1.e-308
MIN_REL_ERROR = 1.e-12

# Helper functions
TWOPI = 2.*math.pi
def mod2pi(f):
    """Returns the angle f modulo 2 pi."""
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def eccentricAnomaly(e,M):
    """Returns the eccentric anomaly given the eccentricity and mean anomaly of a Keplerian orbit.

    Keyword arguments:
    e -- the eccentricity
    M -- the mean anomaly
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
    """Defines the same data structure as in tools.h"""
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
        return "<rebound.Orbit instance, a=%s e=%s>"%(str(self.a),str(self.e))


class Particle(Structure):
    """A particle datastructure. Same as defined in particle.h"""
    _fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double),
                ("vx", c_double),
                ("vy", c_double),
                ("vz", c_double),
                ("ax", c_double),
                ("ay", c_double),
                ("az", c_double),
                ("m", c_double) ]
    def __str__(self):
        return "<rebound.Particle object, m=%s x=%s y=%s z=%s vx=%s vy=%s vz=%s>"%(self.m,self.x,self.y,self.z,self.vx,self.vy,self.vz)
    
    def __init__(self, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None):   
        if particle is not None:
            raise ValueError("Cannot initialise particle from other particles.")
        cart = [x,y,z,vx,vy,vz]
        orbi = [primary,a,anom,e,omega,inc,Omega,MEAN]
        if m is None:   #default value for mass
            m = 0.
        if notNone(cart) and notNone(orbi):
                raise ValueError("You cannot pass cartesian coordinates and orbital elements at the same time.")
        if notNone(orbi):
            if primary is None:
                primary = rebound.module.calculate_com()
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
            self.set_orbit(m=m,primary=primary,a=a,anom=anom,e=e,omega=omega,inc=inc,Omega=Omega,MEAN=MEAN)
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

    def set_orbit(self,
                    m,          # mass
                    primary,    # central body (rebound.Particle object)
                    a,          # semimajor axis
                    anom=0.,    # anomaly
                    e=0.,       # eccentricity
                    omega=0.,   # argument of pericenter
                    inc=0.,     # inclination
                    Omega=0.,   # longitude of ascending node
                    MEAN=False):    # mean anomaly
        """ Initialises a particle structure with the passed set of
            orbital elements. Mass (m), primary and 'a' are required (see Parameters
            below, and any orbital mechanics text, e.g., Murray & Dermott
            Solar System Dynamics for definitions). Other values default to zero.
            All angles should be passed in radians. Units are set by the
            gravitational constant G (default = 1.). If MEAN is set to True, anom is
            taken as the mean anomaly, rather than the true anomaly.
            
            Usage
            _____
            TODO: UPDATE
            primary = rebound.Particle(m=1.) # particle with unit mass at origin & v=0
            
            # test particle (m=0) with specified elements using mean anomaly
            p = kepler_particle(m=0.,primary=primary,a=2.5, anom=math.pi/2,e=0.3,
            omega=math.pi/6,inc=math.pi/3,Omega=0.,MEAN=True)
            
            # m=0.1 particle on circular orbit math.pi/4 from x axis in xy plane
            p = kepler_particle(0.1,primary,2.5,math.pi/4)
            
            Parameters
            __________
            m       : (float)            Mass of the particle
            primary : (rebound.Particle) Particle structure for the central body
            a       : (float)            Semimajor axis
            anom    : (float)            True anomaly (default).
            Mean anomaly if MEAN is set to True
            e       : (float)            Eccentricity
            omega   : (float)            Argument of pericenter
            inc     : (float)            Inclination (to xy plane)
            Omega   : (float)            Longitude of the ascending node
            MEAN    : (boolean)          If False (default), anom = true anomaly
            If True, anom = mean anomaly
            
            Returns
            _______
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
        
        n = math.sqrt(rebound.module.G*(primary.m+m)/(a**3))
        if e>1.:
            v0 = n*a/math.sqrt(-(1.-e**2))
        else:
            v0 = n*a/math.sqrt(1.-e**2)
        
        # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
        self.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
        self.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
        self.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)


    def calculate_orbit(self, primary=None, verbose=False):
        """ Returns a rebound.Orbit object with the keplerian orbital elements
            corresponding to the particle around the central body primary
            (rebound.Particle). Edge cases will return values set to None. If
            verbose is set to True (default=False), error messages are printed
            when a breakout condition is met.
            
            Usage
            _____
            TODO: Update!
            orbit = p2orbit(p,primary)
            print(orbit.e) # gives the eccentricity
            
            orbit = p2orbit(p,primary,verbose=True) # will print out error msgs
            
            Parameters
            __________
            self     : (rebound.Particle) particle for which orbital elements are sought
            primary  : (rebound.Particle) central body
            verbose  : (boolean)          If set to True, will print out error msgs
            
            Returns
            _______
            A rebound.Orbit object (with member variables for the orbital elements)
            """
        if primary is None:
            primary = rebound.module.particles[0]

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
        
        mu = rebound.module.G*(self.m+primary.m)
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

