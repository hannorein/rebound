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

    while rebound.get_t()<0.02*np.pi:
        rebound.step()
    
    Ef = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)


    return [rebound.get_t(), np.fabs((Ef-Ei)/Ei)+1e-16]


N = 200
anoms = np.linspace(-np.pi,np.pi,N)
e0s = np.linspace(0,16,N)

parameters = [(anoms[i], 0.001, e0s[j]) for j in range(N) for i in range(N)]

pool = InterruptiblePool()
res = pool.map(simulation,parameters)
res = np.nan_to_num(res)
megno = np.clip(res[:,0].reshape((N,N)),1.8,4.)
lyaptime = np.absolute(res[:,1].reshape((N,N)))

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

f,axarr = plt.subplots(1)
extent=[anoms.min(), anoms.max(), e0s.min(), e0s.max()]

#for ax in axarr
ax = axarr
ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])
ax.set_xlabel(r"true anomly")
ax.set_ylabel(r"$\xi$  defined as  $e=1-10^{-\xi}$")

#im1 = axarr[0].imshow(megno, vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
#cb1 = plt.colorbar(im1, ax=axarr[0])
#cb1.solids.set_rasterized(True)
#cb1.set_label("MEGNO $\\langle Y \\rangle$")

im2 = axarr.imshow(lyaptime, vmax=1, norm=LogNorm(),aspect='auto', origin="lower", interpolation="nearest", cmap="RdYlGn", extent=extent)
cb2 = plt.colorbar(im2, ax=axarr)
cb2.solids.set_rasterized(True)
cb2.set_label("Energy error")

plt.savefig("2body.pdf")

import os
os.system("open 2body.pdf")
