from ctypes import *
libias15 = CDLL('./libias15.so', RTLD_GLOBAL)

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

def set_G(G):
    c_double.in_dll(libias15, "G").value = G

def set_dt(dt):
    c_double.in_dll(libias15, "dt").value = dt

def set_t(t):
    c_double.in_dll(libias15, "t").value = t

def get_t():
    return c_double.in_dll(libias15, "t").value

def set_particles(particles):
    c_int.in_dll(libias15,"N").value = len(particles)
    arr = (Particle * len(particles))(*particles)
    libias15.setp(byref(arr))

def particles_add(particles):
    if isinstance(particles,list):
        for particle in particles:
            libias15.particles_add(particle)
    else:
       libias15.particles_add(particles)

def move_to_center_of_momentum():
    libias15.tools_move_to_center_of_momentum()

def step():
    libias15.step()

def integrate(tmax):
    libias15.integrate(c_double(tmax))

def get_particle(i):
    N = c_int.in_dll(libias15,"N").value 
    if i>=N:
        return None
    getp = libias15.particle_get
    getp.restype = Particle
    _p = getp(c_int(i))
    return _p

def get_particles():
    N = c_int.in_dll(libias15,"N").value 
    particles = []
    for i in xrange(0,N):
        particles.append(get_particle(i))
    return particles


