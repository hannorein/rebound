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

particles = []
particles.append( Particle(m=1.) )
particles.append( Particle(m=1e-3,x=1.,vy=1.) )

c_double.in_dll(libias15, "G").value = 1
arr = (Particle * len(particles))(*particles)
print POINTER(Particle).in_dll(libias15, "particles")
POINTER(Particle).in_dll(libias15, "particles").value = pointer(arr)
print POINTER(Particle).in_dll(libias15, "particles")
c_int.in_dll(libias15,"N").value = len(particles)

libias15.integrator_ias15_step()

particles = Particle * 2
print particles

print c_double.in_dll(libias15, "G")

