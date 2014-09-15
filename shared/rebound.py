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


def step():
    libias15.integrator_ias15_step()

def get_particles():
    N = c_int.in_dll(libias15,"N").value 
    getp = libias15.getp
    getp.restype = POINTER(Particle)
    _p = getp()
    particles = []
    for i in xrange(0,N):
        particles.append(_p[i])

    return particles


