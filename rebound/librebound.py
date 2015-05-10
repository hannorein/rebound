from ctypes import *
import math
import os
import ctypes.util
import pkg_resources

_pymodulespath = os.path.dirname(__file__)
#Find the rebound C library
try:
    clibrebound = CDLL(_pymodulespath+"/../librebound.so", RTLD_GLOBAL)
except:
    try:
        clibrebound = CDLL(_pymodulespath + '/../shared/librebound/librebound.so', RTLD_GLOBAL)
    except:
        print("Cannot find library 'librebound.so'.")
        raise
    pass

#Make changes for python 2 and 3 compatibility
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x

TINY=1.e-308
MIN_REL_ERROR = 1.e-12


def get_build_str():
    return c_char_p.in_dll(clibrebound, "build_str").value


def status():
    """ Returns a string with a summary of the current status 
        of the simulation
        """
    s= ""
    N = get_N()
    s += "---------------------------------\n"
    s += "Rebound version:     \t" + pkg_resources.require("rebound")[0].version +"\n"
    s += "Build on:            \t" +get_build_str() + "\n"
    s += "Number of particles: \t%d\n" %N       
    s += "Simulation time:     \t%f\n" %get_t()
    if N>0:
        s += "---------------------------------\n"
        p = get_particles()
        for i in range(N):
            s += str(p[i]) + "\n"
    s += "---------------------------------"
    print(s)
# Alias
get_status = status

# Set function pointer for additional forces

AFF = CFUNCTYPE(None)
fp = None
def set_additional_forces(func):
    global fp  # keep references
    fp = AFF(func)
    clibrebound.set_additional_forces(fp)

# Setter/getter of parameters and constants
def set_G(G):
    c_double.in_dll(clibrebound, "G").value = G

def get_G():
    return c_double.in_dll(clibrebound, "G").value

def set_dt(dt):
    c_double.in_dll(clibrebound, "dt").value = dt

def get_dt():
    return c_double.in_dll(clibrebound, "dt").value

def set_t(t):
    c_double.in_dll(clibrebound, "t").value = t

def set_min_dt(t):
    c_double.in_dll(clibrebound, "integrator_ias15_min_dt").value = t

def get_t():
    return c_double.in_dll(clibrebound, "t").value

def init_megno(delta):
    clibrebound.tools_megno_init(c_double(delta))

def get_megno():
    clibrebound.tools_megno.restype = c_double
    return clibrebound.tools_megno()

def get_lyapunov():
    clibrebound.tools_lyapunov.restype = c_double
    return clibrebound.tools_lyapunov()

def get_N():
    return c_int.in_dll(clibrebound,"N").value 

def get_N_megno():
    return c_int.in_dll(clibrebound,"N_megno").value 

def get_iter():
    return c_int.in_dll(clibrebound,"iter").value 

def get_timing():
    return c_double.in_dll(clibrebound,"timing").value 

# Setter functions, particle data
def set_particles(_particles):
    c_int.in_dll(clibrebound,"N").value = len(_particles)
    arr = (Particle * len(_particles))(*_particles)
    clibrebound.setp(byref(arr))

def add_particle(particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None):   
    """Adds a particle to REBOUND. Accepts one of the following four sets of arguments:
    1) A single Particle structure.
    2) A list of Particle structures.
    3) The particle's mass and a set of cartesian coordinates: m,x,y,z,vx,vy,vz.
    3) The primary as a Particle structure, the particle's mass and a set of orbital elements primary,a,anom,e,omega,inv,Omega,MEAN (see kepler_particle() for the definition of orbital elements). 
    """
    if isinstance(particle,Particle):
        clibrebound.particles_add(particle)
    elif isinstance(particle,list):
        for p in particle:
            add_particle(p)
    elif isinstance(particle,str):
        add_particle(horizons.getParticle(**locals()))
    else: 
        add_particle(Particle(**locals()))

# Aliases
add_particles = add_particle
add = add_particle

# Particle getter functions
def get_particle(i):
    N = get_N() 
    if i>=N:
        return None
    getp = clibrebound.particle_get
    getp.restype = Particle
    _p = getp(c_int(i))
    return _p

def get_particles_array():
    N = get_N() 
    _particles = []
    for i in range(0,N):
        _particles.append(particle_get(i))
    return _particles


def get_particles():
    """Return an array that points to the particle structure.
    This is a pointer and thus the contents of the array update 
    as the simulation progresses. Note that the pointer could change,
    for example when a particle is added or removed from the simulation. 
    In that case, get a fresh pointer with get_particles().

    particles() is an alias of this function.
    """
    N = c_int.in_dll(clibrebound,"N").value 
    getp = clibrebound.particles_get
    getp.restype = POINTER(Particle)
    return getp()
# Alias
particles = get_particles

# Orbit getter
def get_orbits(heliocentric=False):
    """ Returns an array of Orbits of length N-1.

    Parameters
    __________

    By default this functions returns the orbits in Jacobi coordinates. 
    Set the parameter heliocentric to True to return orbits in heliocentric coordinates.
    """
    _particles = get_particles()
    orbits = []
    for i in range(1,get_N()):
        if heliocentric:
            com = _particles[0]
        else:
            com = get_com(i)
        orbits.append(_particles[i].get_orbit(primary=com))
    return orbits

# Tools
def move_to_center_of_momentum():
    clibrebound.tools_move_to_center_of_momentum()
# Alias
move_to_barycentric_frame = move_to_center_of_momentum
move_to_center_of_mass = move_to_center_of_momentum
move_to_com = move_to_center_of_momentum

def reset():
    debug.reset_debug()
    clibrebound.reset()

def set_integrator_whfast_corrector(on=11):
    c_int.in_dll(clibrebound, "integrator_whfast_corrector").value = on

def set_integrator(integrator="IAS15"):
    if isinstance(integrator, int):
        clibrebound.integrator_set(c_int(integrator))
        return
    if isinstance(integrator, basestring):
        debug.integrator_fullname = integrator
        debug.integrator_package = "REBOUND"
        if integrator.lower() == "ias15":
            set_integrator(0)
            return
        if integrator.lower() == "whfast":
            set_integrator(1)
            set_integrator_whfast_corrector(11)
            return
        if integrator.lower() == "whfast-nocor":
            set_integrator(1)
            set_integrator_whfast_corrector(0)
            return
        if integrator.lower() == "sei":
            set_integrator(2)
            return
        if integrator.lower() == "wh":
            set_integrator(3)
            return
        if integrator.lower() == "leapfrog":
            set_integrator(4)
            return
        if integrator.lower() == "leap-frog":
            set_integrator(4)
            return
        if integrator.lower() == "hybrid":
            set_integrator(5)
            return
        if integrator.lower() == "mercury":
            debug.integrator_package = "MERCURY"
            return
        if integrator.lower() == "swifter-whm":
            debug.integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-symba":
            debug.integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-helio":
            debug.integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-tu4":
            debug.integrator_package = "SWIFTER"
            return
    raise ValueError("Warning. Intergrator not found.")

def set_integrator_whfast_persistent_particles(is_per=0):
    if isinstance(is_per, int):
        c_int.in_dll(clibrebound, "integrator_whfast_persistent_particles").value = is_per
        return
    raise ValueError("Expecting integer.")

def set_integrator_whfast_synchronize_manually(synchronize_manually=0):
    if isinstance(synchronize_manually, int):
        c_int.in_dll(clibrebound, "integrator_whfast_synchronize_manually").value = synchronize_manually
        return
    raise ValueError("Expecting integer.")

def set_force_is_velocitydependent(force_is_velocitydependent=1):
    if isinstance(force_is_velocitydependent, int):
        c_int.in_dll(clibrebound, "integrator_force_is_velocitydependent").value = force_is_velocitydependent
        return
    raise ValueError("Expecting integer.")

# Input/Output routines

def output_binary(filename):
    clibrebound.output_binary(c_char_p(filename))
save=output_binary
    
def input_binary(filename):
    clibrebound.input_binary(c_char_p(filename))
load=input_binary
    

# Integration
def step():
    clibrebound.rebound_step()

def integrate(tmax,exactFinishTime=1,keepSynchronized=0):
    if debug.integrator_package =="REBOUND":
        clibrebound.integrate(c_double(tmax),c_int(exactFinishTime),c_int(keepSynchronized))
    else:
        debug.integrate_other_package(tmax,exactFinishTime,keepSynchronized)

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


def get_com(last=None):
    """Returns the center of momentum for all particles in the simulation"""
    m = 0.
    x = 0.
    y = 0.
    z = 0.
    vx = 0.
    vy = 0.
    vz = 0.
    ps = get_particles()    # particle pointer
    if last is not None:
        last = min(last,get_N())
    else:
        last = get_N()
    for i in range(last):
    	m  += ps[i].m
    	x  += ps[i].x*ps[i].m
    	y  += ps[i].y*ps[i].m
    	z  += ps[i].z*ps[i].m
    	vx += ps[i].vx*ps[i].m
    	vy += ps[i].vy*ps[i].m
    	vz += ps[i].vz*ps[i].m
    if m>0.:
        x /= m
        y /= m
        z /= m
        vx /= m
        vy /= m
        vz /= m
    return Particle(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)

# Alias
get_center_of_momentum = get_com
get_center_of_mass = get_com

# Import at the end to avoid circular dependence
from .particle import *
from . import horizons
from . import debug
