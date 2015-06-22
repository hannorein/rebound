#!/usr/bin/python
# This example integrates Jupiter and Saturn in the Solar system for a variety of initial conditions.
# Alongside the normal equations of motions, IAS15 is used to integrate the variational equations.
# These can be used to measure the Mean Exponential Growth of Nearby Orbits (MEGNO), a chaos indicator.
# This example script runs 12^2 simulations and plots the MEGNO value. Values close to <Y>=2 correspond 
# to regular quasi-periodic orbits. Higher values of <Y> correspond to chaotic orbits.
 
# Import the rebound module
import rebound
# Import other modules
import numpy as np
import multiprocessing

# Runs one simulation.
def simulation(par):
    saturn_a, saturn_e = par
    rebound.reset()
    rebound.integrator = "whfast"
    rebound.min_dt = 5.
    rebound.dt = 1.
    
    # These parameters are only approximately those of Jupiter and Saturn.
    sun     = rebound.Particle(m=1.)
    rebound.add(sun)
    jupiter = rebound.add(primary=sun,m=0.000954, a=5.204, anom=0.600, omega=0.257, e=0.048)
    saturn  = rebound.add(primary=sun,m=0.000285, a=saturn_a, anom=0.871, omega=1.616, e=saturn_e)

    rebound.move_to_com()
    rebound.init_megno(1e-16)
    rebound.integrate(1e3*2.*np.pi)

    return [rebound.calculate_megno(),1./(rebound.calculate_lyapunov()*2.*np.pi)] # returns MEGNO and Lypunov timescale in years


### Setup grid and run many simulations in parallel
N = 100                      # Grid size, increase this number to see more detail
a = np.linspace(7.,10.,N)   # range of saturn semi-major axis in AU
e = np.linspace(0.,0.5,N)   # range of saturn eccentricity
parameters = []
for _e in e:
    for _a in a:
        parameters.append([_a,_e])


# Run simulations in parallel
pool = rebound.InterruptiblePool()    # Number of threads default to the number of CPUs on the system
print("Running %d simulations on %d threads..." % (len(parameters), pool._processes))
res = np.nan_to_num(np.array(pool.map(simulation,parameters))) 
megno = np.clip(res[:,0].reshape((N,N)),1.8,4.)             # clip arrays to plot saturated 
lyaptimescale = np.clip(np.absolute(res[:,1].reshape((N,N))),1e1,1e5)

### Create plot and save as pdf 
import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Setup plots
f, axarr = plt.subplots(2,figsize=(10,10))
extent = [a.min(), a.max(), e.min(), e.max()]
for ax in axarr:
    ax.set_xlim(extent[0],extent[1])
    ax.set_ylim(extent[2],extent[3])
    ax.set_xlabel("$a_{\mathrm{Saturn}}$ [AU]")
    ax.set_ylabel("$e_{\mathrm{Saturn}}$")


# Plot MEGNO 
im1 = axarr[0].imshow(megno, vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
cb1 = plt.colorbar(im1, ax=axarr[0])
cb1.solids.set_rasterized(True)
cb1.set_label("MEGNO $\\langle Y \\rangle$")

# Plot Lyapunov timescale
im2 = axarr[1].imshow(lyaptimescale, vmin=1e1, vmax=1e5, norm=LogNorm(), aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn", extent=extent)
cb2 = plt.colorbar(im2, ax=axarr[1])
cb2.solids.set_rasterized(True)
cb2.set_label("Lyapunov timescale [years]")

plt.savefig("megno.pdf")

### Automatically open plot (OSX only)
from sys import platform as _platform
if _platform == "darwin":
    import os
    os.system("open megno.pdf")
