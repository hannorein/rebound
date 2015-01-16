from ctypes import *
import numpy
import sys
if sys.version_info < (3,):
    range = xrange

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

def get_t():
    return c_double.in_dll(libias15, "t").value

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

TWOPI = 2.*numpy.pi
def mod2pi(f):
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def get_E(e,M):
    E = M if e<0.8 else numpy.pi
    
    F = E - e*numpy.sin(M) - M
    for i in range(100):
        E = E - F/(1.0-e*numpy.cos(E))
        F = E - e*numpy.sin(E) - M
        if numpy.fabs(F)<1e-16:
            break
    E = mod2pi(E)
    return E

def init_planet(star,       # central object (rebound.Particle object)
                a,          # semimajor axis
                e=None,     # eccentricity
                omega=None, # argument of pericenter
                f=None,     # true anomaly
                inc=None,   # inclination
                Omega=None, # longitude of ascending node
                mass=None,  # planet mass
                M=None):    # mean anomaly
    '''Returns a particle structure to add into rebound for a given set of
        orbital elements. 'star' and 'a' are required, others will default to zero
        (with a warning). Only f OR M should be passed to specify the location of
        the particle on its orbit.'''
    
    # first check for valid inputs and give warnings on defaults
    if e == None: e = 0.; print('e not passed, setting to 0... (circular)')
    if omega == None:
        omega = 0. # only warn if orbit is non-circular
        if e != 0.: print('omega not passed, setting to 0... (pericenter at ascending node)')
    if inc == None: inc = 0.; print('inc not passed, setting to 0... (in xy plane)')
    if Omega == None:
        Omega = 0. # only warn if orbit is inclined
        if inc != 0.: print('Omega not passed, setting to 0... (ascending node along x axis)')
    if f != None and M != None: raise ValueError('Can only pass EITHER the true (f) or the mean (M) anomaly')
    if f == None and M == None:
        f = 0.;
        if e != 0.:
            print('Neither f nor M passed, setting initial angle from pericenter to 0...')
        else:
            print('Neither f nor M passed, setting particle on x axis...')
    if mass == None: mass = 0.; print('mass not passed, setting to 0...')
    
    if not(0. <= e < 1.): raise ValueError('e must be in range [0,1)') # not sure if these equations work for parabolic/hyperbolic obits
    if not(0. <= inc <= numpy.pi): raise ValueError('inc must be in range [0,pi]')
    
    if f == None: # need to calculate f
        E = get_E(e,M)
        f = mod2pi(2.*numpy.arctan(numpy.sqrt((1. + e)/(1. - e))*numpy.tan(0.5*E)))
    
    r = a*(1-e**2)/(1 + e*numpy.cos(f))
    
    # Murray & Dermott Eq. 2.122
    _x  = star.x + r*(numpy.cos(Omega)*numpy.cos(omega+f) - numpy.sin(Omega)*numpy.sin(omega+f)*numpy.cos(inc))
    _y  = star.y + r*(numpy.sin(Omega)*numpy.cos(omega+f) + numpy.cos(Omega)*numpy.sin(omega+f)*numpy.cos(inc))
    _z  = star.z + r*numpy.sin(omega+f)*numpy.sin(inc)
    
    n = numpy.sqrt(get_G()*(star.m+mass)/(a**3))
    
    # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
    _vx = star.vx + (n*a/numpy.sqrt(1-e*e))*((e+numpy.cos(f))*(-numpy.cos(inc)*numpy.cos(omega)*numpy.sin(Omega) - numpy.cos(Omega)*numpy.sin(omega)) - numpy.sin(f)*(numpy.cos(omega)*numpy.cos(Omega) - numpy.cos(inc)*numpy.sin(omega)*numpy.sin(Omega)))
    _vy = star.vy + (n*a/numpy.sqrt(1-e*e))*((e+numpy.cos(f))*(numpy.cos(inc)*numpy.cos(omega)*numpy.cos(Omega) - numpy.sin(omega)*numpy.sin(Omega)) - numpy.sin(f)*(numpy.cos(omega)*numpy.sin(Omega) + numpy.cos(inc)*numpy.cos(Omega)*numpy.sin(omega)))
    _vz = star.vz + (n*a/numpy.sqrt(1-e*e))*((e+numpy.cos(f))*numpy.cos(omega)*numpy.sin(inc) - numpy.sin(f)*numpy.sin(inc)*numpy.sin(omega))
    
    return Particle(m  = mass, x = _x, y = _y, z = _z, vx = _vx, vy = _vy, vz = _vz)

def p2orbit(p, star):
    '''Takes 2 Particle objects (one for object and one for central body,
        which hold Cartesian elements) and returns a Orbit object
        with all the orbital elements'''
    
    if star.m <= 1.e-30: raise ValueError('Star has no mass.')
    o = Orbit()
    mu = get_G()*(p.m+star.m)
    p.x -= star.x
    p.y -= star.y
    p.z -= star.z
    p.vx -= star.vx
    p.vy -= star.vy
    p.vz -= star.vz
    h0 = (p.y*p.vz - p.z*p.vy)                  # angular momentum vector
    h1 = (p.z*p.vx - p.x*p.vz)
    h2 = (p.x*p.vy - p.y*p.vx)
    o.h = numpy.sqrt ( h0*h0 + h1*h1 + h2*h2 )  # abs value of angular moment
    if o.h <= 1.e-30: raise ValueError('Particle orbit is radial.')
    v = numpy.sqrt ( p.vx*p.vx + p.vy*p.vy + p.vz*p.vz )
    o.r = numpy.sqrt ( p.x*p.x + p.y*p.y + p.z*p.z )
    if o.r <= 1.e-30: raise ValueError('Particle and star positions are the same.')
    vr = (p.x*p.vx + p.y*p.vy + p.z*p.vz)/o.r
    e0 = 1./mu*( (v*v-mu/o.r)*p.x - o.r*vr*p.vx )
    e1 = 1./mu*( (v*v-mu/o.r)*p.y - o.r*vr*p.vy )
    e2 = 1./mu*( (v*v-mu/o.r)*p.z - o.r*vr*p.vz )
    o.e = numpy.sqrt( e0*e0 + e1*e1 + e2*e2 )   # eccentricity
    o.a = -mu/( v*v - 2.*mu/o.r )               # semi major axis
    o.P = 2.*numpy.pi*numpy.sqrt( o.a*o.a*o.a/mu )  # period
    o.inc = numpy.arccos( h2/o.h )              # inclination (wrt xy-plane)
    # if pi/2<i<pi it's retrograde
    n0 = -h1                                    # node vector
    n1 =  h0                                    # in xy plane => no z component
    n = numpy.sqrt( n0*n0 + n1*n1 )
    er = p.x*e0 + p.y*e1 + p.z*e2
    if (n<=1.e-30 or o.inc<=1.e-30):            # we are in the xy plane
        o.Omega=0.
        if (o.e <= 1.e-10):                     # omega not defined for circular orbit
            o.omega = 0.
        else:
            if (e1>=0.):
                o.omega=numpy.arccos(e0/o.e)
            else:
                o.omega = 2.*numpy.pi-numpy.arccos(e0/o.e)
    else:
        if (o.e <= 1.e-10):
            o.omega = 0.
        else:
            if (e2>=0.):                        # omega=0 if perictr at asc node
                o.omega=numpy.arccos(( n0*e0 + n1*e1 )/(n*o.e))
            else:
                o.omega=2.*numpy.pi-numpy.arccos(( n0*e0 + n1*e1 )/(n*o.e))
        if (n1>=0.):
            o.Omega = numpy.arccos(n0/n)
        else:
            o.Omega=2.*numpy.pi-numpy.arccos(n0/n)  # Omega=longitude of asc node
    # taken in xy plane from x axis
    
    if (o.e<=1.e-10):                           # circular orbit
        o.f=0.                                  # f has no meaning
        o.l=0.
    else:
        o.f = er/(o.e*o.r)
        ea = (1.-o.r/o.a)/o.e
        if (o.f>1. or o.f<-1.):                 # failsafe
            o.f = numpy.pi - numpy.pi * o.f
            ea  = numpy.pi - numpy.pi * ea
        else:
            o.f = numpy.arccos(o.f)             # true anom=0 if obj at perictr
            ea  = numpy.arccos(ea)              # eccentric anomaly
        
        if (vr<0.):
            o.f=2.*numpy.pi-o.f
            ea =2.*numpy.pi-ea
        
        o.l = ea -o.e*numpy.sin(ea) + o.omega+ o.Omega  # mean longitude
    
    return o


