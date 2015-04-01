# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool


torb = 2.*np.pi
tmax = 200.34476128*torb

def simulation(par):
    anom, dt, e, integrator = par

    e = 1.-pow(10.,-e)
    dt = pow(10.,dt)*torb

    rebound.reset()
    rebound.set_integrator(integrator)
    rebound.set_dt(dt)

    rebound.add_particle(m=1.)
    rebound.add_particle(m=0., x=(1.-e), vy=np.sqrt((1.+e)/(1.-e)))
    particles = rebound.get_particles()
    
    rebound.init_integrator()
    
    Ei = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)

    rebound.integrate(tmax,0)
    
    Ef = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)

    return [float(rebound.get_iter())/rebound.get_t()*dt, np.fabs((Ef-Ei)/Ei)+1e-16, rebound.get_timing()/rebound.get_t()*dt]


N = 20
dts = np.linspace(-4,-1,N)
e0s = np.linspace(0,6,N)
integrators= ["wh","mikkola"]

niter = []
energyerror = []
timing = []

for integrator in integrators:
    parameters = [(0., dts[i], e0s[j], integrator) for j in range(N) for i in range(N)]

    pool = InterruptiblePool(16)
    res = pool.map(simulation,parameters)
    res = np.nan_to_num(res)
    niter.append(res[:,0].reshape((N,N)))
    energyerror.append(res[:,1].reshape((N,N)))
    timing.append(res[:,2].reshape((N,N)))

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm

f,axarr = plt.subplots(3,2,figsize=(15,15))
extent=[dts.min(), dts.max(), e0s.min(), e0s.max()]

for ay in axarr:
    for ax in ay:
        # ax = axarr
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_xlabel(r"log10$(dt/t_{orb})$")
        ax.set_ylabel(r"$\xi$, where $e=1-10^{-\xi}$")

for i, integrator in enumerate(integrators):
    im1 = axarr[0,i].imshow(energyerror[i], norm=LogNorm(), vmax=np.max(energyerror), vmin=1e-16, aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
    cb1 = plt.colorbar(im1, ax=axarr[0,i])
    cb1.solids.set_rasterized(True)
    cb1.set_label("Relative energy error, " +integrator)

    im2 = axarr[1,i].imshow(niter[i], vmin=0, vmax=np.max(niter), aspect='auto', origin="lower", interpolation="nearest", cmap="RdYlGn_r", extent=extent)
    cb2 = plt.colorbar(im2, ax=axarr[1,i])
    cb2.solids.set_rasterized(True)
    cb2.set_label("Number of iterations (neg = bisection), " +integrator)

    im3 = axarr[2,i].imshow(timing[i], norm=LogNorm(), vmin=np.mean(timing)/4., vmax=4.*np.median(timing), aspect='auto', origin="lower", interpolation="nearest", cmap="RdYlGn_r", extent=extent)
    cb3 = plt.colorbar(im3, ax=axarr[2,i])
    cb3.solids.set_rasterized(True)
    cb3.set_label("Runtime per step (s), " +integrator)
    tick_locator = ticker.LogLocator(base=2.)
    cb3.locator = tick_locator
    cb3.update_ticks()
    plt.show()


plt.savefig("2body.pdf")
print "Average speedup (WH/Mikkola): %.4f" %(np.mean(timing[0])/np.mean(timing[1]))
import os
os.system("open 2body.pdf")
