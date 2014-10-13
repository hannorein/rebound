from ctypes import *
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
            print "Cannot find library 'libias15.so'. Check path set in 'rebound.py'."
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

# Defines the same datastructure as in orbit.h
class Orbit(Structure):
    _fields_ = [("a", c_double),
                ("r", c_double), # Radial distance from central object
                ("h", c_double), # Angular momentum
                ("P", c_double), # Orbital period
                ("l", c_double),
                ("e", c_double),
                ("inc", c_double),
                ("Omega", c_double), # longitude of ascending node
                ("omega", c_double), # argument of perihelion
                ("f", c_double) ]    # true anomaly

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

def set_dt(dt):
    c_double.in_dll(libias15, "dt").value = dt

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

def particle_add(particles):
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
    for i in xrange(0,N):
        particles.append(particle_get(i))
    return particles


def particles_get():
    N = c_int.in_dll(libias15,"N").value
    getp = libias15.particles_get
    getp.restype = POINTER(Particle)
    return getp()


# Tools
def tools_p2orbit(p, star):
    tools_p2orbit = libias15.tools_p2orbit
    tools_p2orbit.restype = Orbit
    orbit = libias15.tools_p2orbit(p, star)
    return orbit

def tools_get_center_of_mass(p1, p2):
    tools_get_center_of_mass = libias15.tools_get_center_of_mass
    tools_get_center_of_mass.restype = Particle
    particle = libias15.tools_get_center_of_mass(p1, p2)
    return particle

def move_to_center_of_momentum():
    libias15.tools_move_to_center_of_momentum()

# Integration
def step():
    libias15.step()

def integrate(tmax):
    libias15.integrate(c_double(tmax))


