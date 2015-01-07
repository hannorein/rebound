#!/usr/bin/python
# This example integrates Jupiter and Saturn in the Solar system for a variety of initial conditions.
# Alongside the normal equations of motions, IAS15 is used to integrate the variational equations.
# These can be used to measure the Mean Exponential Growth of Nearby Orbits (MEGNO), a chaos indicator.
# This example script runs 12^2 simulations and plots the MEGNO value. Values close to <Y>=2 correspond 
# to regular quasi-periodic orbits. Higher values of <Y> correspond to chaotic orbits.
 
# Import the rebound module
import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import math
import os
from interruptible_pool import InterruptiblePool

# mod2pi helper function.
TWOPI = 2.*math.pi
def mod2pi(f):
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

# Convert from orbital parameters to cartesian coordinates
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
    return rebound.Particle(
        m  = mass, 
        x  = cosomega*_x  - sinomega*_y, y  = sinomega*_x  + cosomega*_y, z  = 0,
        vx = cosomega*_vx - sinomega*_vy, vy = sinomega*_vx + cosomega*_vy, vz = 0.)

# Runs one simulation.
def megno(par):
    saturn_P, saturn_e = par
    rebound.reset()
    rebound.set_min_dt(0.1)
    
    # These parameters are only approximately those od Jupiter and Saturn.
    sun     = rebound.Particle( m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    jupiter = setup_planet(sun,0.00095479194, 74.5366, 0.600331, 0.2570604, 0.04838624)
    saturn  = setup_planet(sun,0.00028588598, saturn_P*TWOPI, 0.871866, 1.6161553, saturn_e)

    rebound.particle_add(sun)
    rebound.particle_add(jupiter) 
    rebound.particle_add(saturn) 
    rebound.move_to_center_of_momentum()
    rebound.megno_init(1e-16)

    rebound.integrate(1e3*2.*math.pi)

    return rebound.get_megno()


### Setup grid and run many simulations in parallel
N = 12                      # Grid size, increase this number to see more detail
P = np.linspace(23.,31.,N)  # range of saturn period in years
e = np.linspace(0.,0.5,N)   # range of saturn eccentricity
v = []
for _e in e:
    for _P in P:
        v.append([_P,_e])

pool = InterruptiblePool(24)    # Number of threads
res = pool.map(megno,v)         # Run simulations in parallel

### Create plot and save as pdf 
import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as pl

extent = [P.min(), P.max(), e.min(), e.max()]

pl.imshow(np.array(res).reshape((N,N)), vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
cb1 = pl.colorbar()
cb1.solids.set_rasterized(True)
cb1.set_label("MEGNO stability indicator $\\langle Y \\rangle$")
pl.xlim(extent[0],extent[1])
pl.ylim(extent[2],extent[3])
pl.xlabel("$P_{\mathrm{Saturn}}$ [yrs]")
pl.ylabel("$e_{\mathrm{Saturn}}$")

pl.savefig("megno.pdf")

### Show plot (OSX only)
os.system("open megno.pdf")
