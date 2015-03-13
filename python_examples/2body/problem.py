# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool

def simulation(par):
    anom, dt, e = par

    e = 1.-pow(10.,-e);
    rebound.reset()
    rebound.set_integrator("mikkola")
    rebound.set_dt(dt)

    rebound.particle_add(m=1.)
    rebound.particle_add(m=0., a=1., e=e, anom=anom)
    particles = rebound.particles_get()

    Ei = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)

    #while rebound.get_t()<0.2*np.pi:
    rebound.step()
    
    Ef = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)


    return [rebound.get_iter(), np.fabs((Ef-Ei)/Ei)+1e-16]


N = 200
anoms = np.linspace(-np.pi,np.pi,N)
e0s = np.linspace(0,14,N)

parameters = [(anoms[i], 0.01*2.*np.pi, e0s[j]) for j in range(N) for i in range(N)]

pool = InterruptiblePool()
res = pool.map(simulation,parameters)
res = np.nan_to_num(res)
niter = res[:,0].reshape((N,N))
energyerror = res[:,1].reshape((N,N))

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

f,axarr = plt.subplots(2,figsize=(10,10))
extent=[anoms.min(), anoms.max(), e0s.min(), e0s.max()]

for ax in axarr:
    # ax = axarr
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel(r"true anomly")
    ax.set_ylabel(r"$\xi$  defined as  $e=1-10^{-\xi}$")

im1 = axarr[0].imshow(energyerror, norm=LogNorm(), aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
cb1 = plt.colorbar(im1, ax=axarr[0])
cb1.solids.set_rasterized(True)
cb1.set_label("Relative energy error")

im2 = axarr[1].imshow(niter, vmin=-3, aspect='auto', origin="lower", interpolation="nearest", cmap="RdYlGn", extent=extent)
cb2 = plt.colorbar(im2, ax=axarr[1])
cb2.solids.set_rasterized(True)
cb2.set_label("Number of iterations (neg = bisection)")

plt.savefig("2body.pdf")

import os
os.system("open 2body.pdf")
