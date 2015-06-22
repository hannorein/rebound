#!/usr/bin/python
import rebound
from rebound.interruptible_pool import InterruptiblePool
# Import other modules
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.colors import LogNorm

# Runs one simulation.
def simulation(par):
    saturn_a, saturn_e = par
    rebound.reset()
    rebound.integrator = "whfast"
    rebound.integrator_whfast_safe_mode = 0 
    rebound.dt = 5.
    
    # These parameters are only approximately those of Jupiter and Saturn.
    rebound.add(m=1.)
    rebound.add(m=0.000954, a=5.204, anom=0.600, omega=0.257, e=0.048)
    rebound.add(m=0.000285, a=saturn_a, anom=0.871, omega=1.616, e=saturn_e)

    rebound.move_to_com()
    rebound.init_megno(1e-16)
    rebound.integrate(5e2*2.*np.pi) # integrator for 500 years

    return [rebound.calculate_megno(),1./(rebound.calculate_lyapunov()*2.*np.pi)] # returns MEGNO and Lypunov timescale in years


def updatePlot(first=False):
    # This constructs a 2d array.
    # The current implementation is slow, but simple.
    keys = np.array(resd.keys())
    x1 = np.unique(keys.T[0])
    x2 = np.unique(keys.T[1])
    res = np.empty((len(x2),len(x1),2))
    for i,_x1 in enumerate(x1):
        for j,_x2 in enumerate(x2):
            res[j][i] = resd[(_x1,_x2)]

    # Clip arrays 
    megno = np.clip(res[:,:,0],1.8,4.)             
    lyaptimescale = np.clip(np.absolute(res[:,:,1]),1e1,4e3)

    # Plot MEGNO 
    im1 = axarr[0].imshow(megno, vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)

    # Plot Lyapunov timescale
    im2 = axarr[1].imshow(lyaptimescale, vmin=1e1, vmax=4e3, norm=LogNorm(), aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn", extent=extent)

    if first:
        cb1 = plt.colorbar(im1, ax=axarr[0])
        cb1.solids.set_rasterized(True)
        cb1.set_label("MEGNO $\\langle Y \\rangle$")

        cb2 = plt.colorbar(im2, ax=axarr[1])
        cb2.solids.set_rasterized(True)
        cb2.set_label("Lyapunov timescale [years]")
    plt.draw()

pool = InterruptiblePool()    # Number of threads default to the number of CPUs on the system
def runSim(p):
    print("Running %d simulations." % len(p))
    res = np.nan_to_num(np.array(pool.map(simulation,p)))
    for i,r in enumerate(res):
        resd[p[i]] = r 

# Setup grid and run many simulations in parallel
a = np.array([7.,10.])   # range of saturn semi-major axis in AU
e = np.array([0.,0.5])   # range of saturn eccentricity

# Setup plots
f, axarr = plt.subplots(2,figsize=(10,8))
extent = [a.min(), a.max(), e.min(), e.max()]
for ax in axarr:
    ax.set_xlim(extent[0],extent[1])
    ax.set_ylim(extent[2],extent[3])
    ax.set_xlabel("$a_{\mathrm{Saturn}}$ [AU]")
    ax.set_ylabel("$e_{\mathrm{Saturn}}$")

# Results are stored in this dictionary
resd = {}

# Initial parameters (2x2 grid)
parameters = [(_a, _e) for _a in a for _e in e]

# Run and plot first simulations
runSim(parameters)
updatePlot(first=True)

# Eight levels of refinement
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
