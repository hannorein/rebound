from ctypes import *
import math
try:
    range = xrange  # this means we're using python 2.x
except NameError:
    pass            # this means we're using python >2.x

# Try to load libias15 from the obvioud places it could be in.
try:
    libias15 = CDLL('../../shared/libias15.so', RTLD_GLOBAL)
except:
    try:
        libias15 = CDLL('../shared/libias15.so', RTLD_GLOBAL)
    except:
        try:
            libias15 = CDLL('shared/libias15.so', RTLD_GLOBAL)
        except:
            print("Cannot find library 'libias15.so'. Check path set in 'rebound.py'.")
            raise


# Defines the same datastructure as in particle.h
class Particle(Structure):
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

# Defines the same data structure as in tools.h
class Orbit():
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

# Set function pointer for additional forces

AFF = CFUNCTYPE(None)
fp = None
def set_additional_forces(func):
    global fp  # keep references
    fp = AFF(func)
    libias15.set_additional_forces(fp)

# Setter/getter of parameters and constants
def set_G(G):
    c_double.in_dll(libias15, "G").value = G

def get_G():
    return c_double.in_dll(libias15, "G").value

def set_dt(dt):
    c_double.in_dll(libias15, "dt").value = dt

def get_dt():
    return c_double.in_dll(libias15, "dt").value

def set_t(t):
    c_double.in_dll(libias15, "t").value = t

def set_min_dt(t):
    c_double.in_dll(libias15, "integrator_min_dt").value = t

def get_t():
    return c_double.in_dll(libias15, "t").value

def megno_init(delta):
    libias15.integrator_megno_init(c_double(delta))

def get_megno():
    libias15.integrator_megno.restype = c_double
    return libias15.integrator_megno()

def get_N():
    return c_int.in_dll(libias15,"N").value 


# Setter/getter of particle data
def set_particles(particles):
    c_int.in_dll(libias15,"N").value = len(particles)
    arr = (Particle * len(particles))(*particles)
    libias15.setp(byref(arr))

def particles_add(particles):
    particle_add(particles)

def particle_add(particles=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None):
    if any([m,x,y,z,vx,vy,vz]):
        if particles is not None:
            raise ValueError("You cannot have a particle structure and float values at the same time.")
        if m is None:
            m = 0.
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
        particles = Particle(m=m,x=x,y=y,z=z,vx=vx,vy=vy,vz=vz)
    if isinstance(particles,list):
        for particle in particles:
            libias15.particles_add(particle)
    else:
       libias15.particles_add(particles)

def particle_get(i):
    N = get_N() 
    if i>=N:
        return None
    getp = libias15.particle_get
    getp.restype = Particle
    _p = getp(c_int(i))
    return _p

def particles_get_array():
    N = get_N() 
    particles = []
    for i in range(0,N):
        particles.append(particle_get(i))
    return particles


def particles_get():
    N = c_int.in_dll(libias15,"N").value 
    getp = libias15.particles_get
    getp.restype = POINTER(Particle)
    return getp()


# Tools
def move_to_center_of_momentum():
    libias15.tools_move_to_center_of_momentum()

def reset():
    libias15.reset()

# Integration
def step():
    libias15.ias15_step()

def integrate(tmax):
    libias15.integrate(c_double(tmax))

TWOPI = 2.*math.pi
def mod2pi(f):
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def get_E(e,M):
    E = M if e<0.8 else math.pi
    
    F = E - e*math.sin(M) - M
    for i in range(100):
        E = E - F/(1.0-e*math.cos(E))
        F = E - e*math.sin(E) - M
        if math.fabs(F)<1e-16:
            break
    E = mod2pi(E)
    return E

def kepler_particle(m,
                    primary,    # central body (rebound.Particle object)
                    a,          # semimajor axis
                    anom=0.,  # anomaly
                    e=0.,     # eccentricity
                    omega=0., # argument of pericenter
                    inc=0.,   # inclination
                    Omega=0., # longitude of ascending node
                    MEAN=False):    # mean anomaly
    """Returns a particle structure initialized with the passed set of
        orbital elements. Mass (m), primary and 'a' are required (see Parameters
        below, and any orbital mechanics text, e.g., Murray & Dermott
        Solar System Dynamics for definitions). Other values default to zero.
        All angles should be passed in radians. Units are set by the
        gravitational constant G (default = 1.). If MEAN is set to True, anom is
        taken as the mean anomaly, rather than the true anomaly.
        
        Usage
        _____
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
    
    if not(0.<=e<1.): raise ValueError('e must be in range [0,1)')
    # not sure if these equations work for parabolic/hyperbolic obits
    if not(0.<=inc<=math.pi): raise ValueError('inc must be in range [0,pi]')
    
    if MEAN is True: # need to calculate f
        E = get_E(e,anom)
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
    x  = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
    y  = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
    z  = primary.z + r*(so*cf+co*sf)*si
    
    n = math.sqrt(get_G()*(primary.m+m)/(a**3))
    v0 = n*a/math.sqrt(1.-e**2)
    
    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
    vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
    vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)
    
    return Particle(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

TINY=1.e-308
MIN_REL_ERROR = 1.e-12
# from tools.c.  Converts cartesian elements to orbital elements.

def p2orbit(p, primary,verbose=False):
    """ Returns a rebound.Orbit object with the keplerian orbital elements
        corresponding to p (rebound.Particle) around the central body primary
        (rebound.Particle). Edge cases will return values set to None. If
        verbose is set to True (default=False), error messages are printed
        when a breakout condition is met.
        
        Usage
        _____
        orbit = p2orbit(p,primary)
        print(orbit.e) # gives the eccentricity
        
        orbit = p2orbit(p,primary,verbose=True) # will print out error msgs
        
        Parameters
        __________
        p        : (rebound.Particle) particle for which orbital elements are sought
        primary  : (rebound.Particle) central body
        verbose  : (boolean)          If set to True, will print out error msgs
        
        Returns
        _______
        A rebound.Orbit object (with member variables for the orbital elements)
        """
    o = Orbit()
    if primary.m <= TINY:
        if verbose is True:
            print("Star has no mass.")
        return o                            # all values set to None
    
    dx = p.x - primary.x
    dy = p.y - primary.y
    dz = p.z - primary.z
    o.r = math.sqrt ( dx*dx + dy*dy + dz*dz )
    if o.r <= TINY:
        if verbose is True:
            print('Particle and primary positions are the same.')
        return o
    
    dvx = p.vx - primary.vx
    dvy = p.vy - primary.vy
    dvz = p.vz - primary.vz
    v = math.sqrt ( dvx*dvx + dvy*dvy + dvz*dvz )
    
    mu = get_G()*(p.m+primary.m)
    o.a = -mu/( v*v - 2.*mu/o.r )               # semi major axis
    
    h0 = (dy*dvz - dz*dvy)                      # angular momentum vector
    h1 = (dz*dvx - dx*dvz)
    h2 = (dx*dvy - dy*dvx)
    o.h = math.sqrt ( h0*h0 + h1*h1 + h2*h2 )   # abs value of angular momentum
    if o.h/o.r/v <= MIN_REL_ERROR:
        if verbose is True:
            print('Particle orbit is radial.')
        return o
    
    vr = (dx*dvx + dy*dvy + dz*dvz)/o.r
    e0 = 1./mu*( (v*v-mu/o.r)*dx - o.r*vr*dvx )
    e1 = 1./mu*( (v*v-mu/o.r)*dy - o.r*vr*dvy )
    e2 = 1./mu*( (v*v-mu/o.r)*dz - o.r*vr*dvz )
    o.e = math.sqrt( e0*e0 + e1*e1 + e2*e2 )   # eccentricity
    
    o.P = 2.*math.pi*math.sqrt( o.a*o.a*o.a/mu )  # period
    o.inc = math.acos( h2/o.h )               # inclination (wrt xy-plane)
    # if pi/2<i<pi it's retrograde
    n0 = -h1                                    # node vector
    n1 =  h0                                    # in xy plane => no z component
    n = math.sqrt( n0*n0 + n1*n1 )
    er = dx*e0 + dy*e1 + dz*e2
    if (n/o.r/v<=MIN_REL_ERROR or o.inc<=MIN_REL_ERROR):# we are in the xy plane
        o.Omega=0.
        if (o.e <= MIN_REL_ERROR):              # omega not defined for circular orbit
            o.omega = 0.
        else:
            if (e1>=0.):
                o.omega=math.acos(e0/o.e)
            else:
                o.omega = 2.*math.pi-math.acos(e0/o.e)
    else:
        if (o.e <= MIN_REL_ERROR):
            o.omega = 0.
        else:
            if (e2>=0.):                        # omega=0 if perictr at asc node
                o.omega=math.acos(( n0*e0 + n1*e1 )/(n*o.e))
            else:
                o.omega=2.*math.pi-math.acos(( n0*e0 + n1*e1 )/(n*o.e))
        if (n1>=0.):
            o.Omega = math.acos(n0/n)
        else:
            o.Omega=2.*math.pi-math.acos(n0/n)# Omega=longitude of asc node
    # taken in xy plane from x axis
    
    if (o.e<=MIN_REL_ERROR):                           # circular orbit
        o.f=0.                                  # f has no meaning
        o.l=0.
    else:
        o.f = er/(o.e*o.r)
        ea = (1.-o.r/o.a)/o.e
        if (o.f>1. or o.f<-1.):                 # failsafe
            o.f = math.pi/2. - math.pi/2.*o.f
            ea  = math.pi/2. - math.pi/2.*ea
        else:
            o.f = math.acos(o.f)             # true anom=0 if obj at perictr
            ea  = math.acos(ea)              # eccentric anomaly
        
        if (vr<0.):
            o.f=2.*math.pi-o.f
            ea =2.*math.pi-ea
        
        o.l = ea -o.e*math.sin(ea) + o.omega+ o.Omega  # mean longitude
    
    return o

