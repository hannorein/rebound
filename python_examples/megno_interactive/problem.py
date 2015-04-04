#!/usr/bin/python
# Import the rebound module
import sys; sys.path.append('../../python_modules')
import rebound
from interruptible_pool import InterruptiblePool
# Import other modules
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.colors import LogNorm

# Runs one simulation.
def simulation(par):
    saturn_a, saturn_e = par
    #/return [saturn_a, saturn_e]
    rebound.reset()
    rebound.set_integrator("mikkola")
    rebound.set_masses_are_constant(1)
    rebound.set_min_dt(5.)
    rebound.set_dt(1.)
    
    # These parameters are only approximately those of Jupiter and Saturn.
    sun     = rebound.Particle(m=1.)
    rebound.add_particle(sun)
    jupiter = rebound.add_particle(primary=sun,m=0.000954, a=5.204, anom=0.600, omega=0.257, e=0.048)
    saturn  = rebound.add_particle(primary=sun,m=0.000285, a=saturn_a, anom=0.871, omega=1.616, e=saturn_e)

    rebound.move_to_center_of_momentum()
    rebound.init_megno(1e-16)
    rebound.integrate(1e3*2.*np.pi)

    return [rebound.get_megno(),1./(rebound.get_lyapunov()*2.*np.pi)] # returns MEGNO and Lypunov timescale in years


pool = InterruptiblePool()    # Number of threads default to the number of CPUs on the system

### Setup grid and run many simulations in parallel
a = np.array([7.,10.])   # range of saturn semi-major axis in AU
e = np.array([0.,0.5])   # range of saturn eccentricity

# Setup plots
f, axarr = plt.subplots(2,figsize=(10,10))
extent = [a.min(), a.max(), e.min(), e.max()]
for ax in axarr:
    ax.set_xlim(extent[0],extent[1])
    ax.set_ylim(extent[2],extent[3])
    ax.set_xlabel("$a_{\mathrm{Saturn}}$ [AU]")
    ax.set_ylabel("$e_{\mathrm{Saturn}}$")

#cb1 = plt.colorbar(im1, ax=axarr[0])
#cb1.solids.set_rasterized(True)
#cb1.set_label("MEGNO $\\langle Y \\rangle$")
#
#cb2 = plt.colorbar(im2, ax=axarr[1])
#cb2.solids.set_rasterized(True)
#cb2.set_label("Lyapunov timescale [years]")

parameters = [(_a, _e) for _a in a for _e in e]
resd = {}
res = np.nan_to_num(np.array(pool.map(simulation,parameters)))
for i,r in enumerate(res):
    resd[parameters[i]] = r 

def toArray(d):
    keys = np.array(d.keys())
    x1 = np.unique(keys.T[0])
    x2 = np.unique(keys.T[1])
    c = np.empty((len(x2),len(x1),2))
    for i,_x1 in enumerate(x1):
        for j,_x2 in enumerate(x2):
            c[j][i] = d[(_x1,_x2)]
    return c

def updatePlot():
    res = toArray(resd)
    # Run simulations in parallel
    megno = np.clip(res[:,:,0],1.8,4.)             # clip arrays to plot saturated 
    lyaptimescale = np.clip(np.absolute(res[:,:,1]),1e1,1e5)

    # Plot MEGNO 
    im1 = axarr[0].imshow(megno, vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)

    # Plot Lyapunov timescale
    im2 = axarr[1].imshow(lyaptimescale, vmin=1e1, vmax=1e5, norm=LogNorm(), aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn", extent=extent)
    plt.draw()


def runSim(p):
    print("Running %d simulations." % len(p))
    res = np.nan_to_num(np.array(pool.map(simulation,p)))
    for i,r in enumerate(res):
        resd[p[i]] = r 

for i in xrange(8):
    _a = np.linspace((a[0]+a[1])/2.,a[-1],len(a))[:-1]
    a = np.sort(np.concatenate((a,_a)))
    parameters = [(__a, _e) for __a in _a for _e in e]
    runSim(parameters)
    updatePlot()
    
    _e = np.linspace((e[0]+e[1])/2.,e[-1],len(e))[:-1]
    e = np.sort(np.concatenate((e,_e)))
    parameters = [(_a, __e) for _a in a for __e in _e]
    runSim(parameters)
    updatePlot()


raw_input('Press enter...')
