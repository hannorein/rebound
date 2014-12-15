# Import the rebound module
import sys; sys.path.append('../')
import rebound
import numpy as np
from rebound import Particle
import math
from interruptible_pool import InterruptiblePool

TWOPI = 2.*math.pi
def mod2pi(f):
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def setup_planet(com, mass, period, M, omega, eccentricity):
    mu = com.m+mass
    n = TWOPI/period
    a = pow(mu/(n*n),1./3.)
    E = M
    if eccentricity>0.8:
        E = np.pi
    F = E - eccentricity*math.sin(M) - M
    for i in xrange(100):
        E = E - F/(1.0-eccentricity*math.cos(E))
        F = E - eccentricity*math.sin(E) - M
        if math.fabs(F)<1e-16:
            break
    E = mod2pi(E)
    f = mod2pi(2.*math.atan(math.sqrt((1. + eccentricity)/(1. - eccentricity))*math.tan(0.5*E)))
    r = a*(1.-eccentricity*eccentricity)/(1.+eccentricity*math.cos(f))
    n = TWOPI/period
    _x  = com.x  + r* math.cos(f)
    _y  = com.y  + r* math.sin(f)
    _vx = com.vx - n*a/math.sqrt(1.-eccentricity*eccentricity)*math.sin(f)
    _vy = com.vy + n*a/math.sqrt(1.-eccentricity*eccentricity)*(eccentricity+math.cos(f))
    cosomega = math.cos(omega)
    sinomega = math.sin(omega)
    return Particle(
        m  = mass, 
        x  = cosomega*_x  - sinomega*_y,
        y  = sinomega*_x  + cosomega*_y,
        z  = 0,
        vx = cosomega*_vx - sinomega*_vy,
        vy = sinomega*_vx + cosomega*_vy,
        vz = 0.)

def megno(par):
    saturn_P, saturn_e = par
    rebound.reset()
    rebound.set_min_dt(0.1)
    sun = Particle( m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    jupiter = setup_planet(sun,0.00095479194, 74.5366, 0.600331, 0.2570604, 0.04838624)
    saturn  = setup_planet(sun,0.00028588598, saturn_P, 0.871866, 1.6161553, saturn_e)

    rebound.particle_add(sun)
    rebound.particle_add(jupiter) 
    rebound.particle_add(saturn) 
    rebound.move_to_center_of_momentum()
    rebound.megno_init()

    rebound.integrate(1e3*2.*math.pi)

    return rebound.get_megno()

N = 128
grid = np.empty([N*N,2])
v = []
a = np.linspace(150.,200.,N)
e = np.linspace(0.,0.5,N)
for _a in a:
    for _e in e:
        v.append([_a,_e])


pool = InterruptiblePool(4)  # 2 threads
res = pool.map(megno,v)

heatmap = np.zeros((N,N,3))
for i, _a in enumerate(a):
    for j, _e in enumerate(e):
        heatmap[i][j][0] = _e
        heatmap[i][j][1] = _a
        heatmap[i][j][2] = res[j*N+i]
np.save("megno.npy",heatmap)
