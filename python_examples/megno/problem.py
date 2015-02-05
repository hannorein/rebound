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
import os
import multiprocessing

# Runs one simulation.
def megno(par):
    saturn_a, saturn_e = par
    rebound.reset()
    rebound.set_min_dt(0.1)
    
    # These parameters are only approximately those od Jupiter and Saturn.
    sun     = rebound.Particle(m=1.)
    rebound.particle_add(sun)
    jupiter = rebound.particle_add(primary=sun,m=0.000954, a=5.204, anom=0.600, omega=0.257, e=0.048)
    saturn  = rebound.particle_add(primary=sun,m=0.000285, a=saturn_a, anom=0.871, omega=1.616, e=saturn_e)

    rebound.move_to_center_of_momentum()
    rebound.megno_init(1e-16)

    rebound.integrate(1e3*2.*np.pi)

    return rebound.get_megno()


### Setup grid and run many simulations in parallel
N = 12                      # Grid size, increase this number to see more detail
a = np.linspace(7.,10.,N)   # range of saturn semi-major axis in AU
e = np.linspace(0.,0.5,N)   # range of saturn eccentricity
v = []
for _e in e:
    for _a in a:
        v.append([_a,_e])

pool = multiprocessing.Pool(24)    # Number of threads
res = pool.map(megno,v)            # Run simulations in parallel

### Create plot and save as pdf 
import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as pl

extent = [a.min(), a.max(), e.min(), e.max()]

pl.imshow(np.array(res).reshape((N,N)), vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
cb1 = pl.colorbar()
cb1.solids.set_rasterized(True)
cb1.set_label("MEGNO stability indicator $\\langle Y \\rangle$")
pl.xlim(extent[0],extent[1])
pl.ylim(extent[2],extent[3])
pl.xlabel("$a_{\mathrm{Saturn}}$ [AU]")
pl.ylabel("$e_{\mathrm{Saturn}}$")

pl.savefig("megno.pdf")

### Show plot (OSX only)
os.system("open megno.pdf")
