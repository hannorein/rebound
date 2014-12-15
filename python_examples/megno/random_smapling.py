# Import the rebound module
import sys; sys.path.append('../')
import rebound
import numpy as np
from rebound import Particle
import math
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

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

pool = InterruptiblePool(4)  # 2 threads
plt.ion()
plt.figure()

N = 128
xmin = 150.
xmax = 200.
ymin = 0.0
ymax = 0.5
def normalize_x(data):
    data = data.astype(np.float)
    return (data - xmin) / (xmax - xmin)

def normalize_y(data):
    data = data.astype(np.float)
    return (data - ymin) / (ymax - ymin)

xi = np.linspace(xmin,xmax,N)
yi = np.linspace(ymin,ymax,N)
xi, yi = np.meshgrid(xi, yi)
    
xi_new = normalize_x(xi)
yi_new = normalize_y(yi)

x = np.array([])
y = np.array([])
z = np.array([])

first = 1
for j in xrange(100):
    v = []
    for i in xrange(16): # random samples
        a = np.random.uniform(xmin,xmax)
        e = np.random.uniform(ymin,ymax)
        v.append([a,e])

    res = pool.map(megno,v)

    x = np.hstack((x,np.array(v)[:,0]))
    y = np.hstack((y,np.array(v)[:,1]))
    z = np.hstack((z,np.array(res)))


    x_new = normalize_x(x)
    y_new = normalize_y(y)
    zi = mlab.griddata(x_new, y_new, z, xi_new, yi_new)

    plt.pcolormesh(xi,yi,zi,vmin=0,vmax=4.,cmap="Blues")
    if first:
        first = 0
        plt.colorbar()
        plt.axis([xmin, xmax, ymin, ymax])
    plt.draw()
plt.show()


plt.show()

