# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool

def simulation(par):
	S, dt,e0 = par
	
	rebound.reset()
	rebound.set_integrator("mikkola")
	rebound.set_dt(dt)

	r0 = 1.-e0 # assumes the orbit has oculating a = 1 and we start at pericenter
	vy = np.sqrt(2./(1.-e0)-1.)
	
	rebound.particle_add(m=1.)
	rebound.particle_add(m=0.,a=1.,e=e0,omega=0.,inc=0.,Omega=0.,MEAN=True)
	
	rebound.move_to_center_of_momentum()
	rebound.megno_init(1.e-16)
	
	particles = rebound.particles_get()

	def starkforce(): # need to put inside simulation(par) to have access to S and particles
		particles[1].ax += S
	
	rebound.set_additional_forces(starkforce)
	
	rebound.integrate(50.*np.pi)

	'''	
	xs = []
	ys = []	
	steps = 0
	while rebound.get_t()<14*np.pi:
		rebound.step()
		steps += 1
		if steps % 1 == 0:
			xs.append(particles[1].x)
			ys.append(particles[1].y)
	
	n=4000
	fig, ax = plt.subplots()
	ax.plot(xs[::],ys[::])
	plt.show()
	'''
	
	return [rebound.get_megno(), rebound.get_t()]

#I always set the (osculating) semimajor axis to 1, you can pass different initial e values

e0 = 0.9 # Rauch uses 0.9 for Fig 4
Scrit = 0.25 # always true if you use G=M=a=1

N = 200
dts = np.linspace(0.01,1.,N)
Ss = np.linspace(0,0.5,N)
parameters = [(Ss[i]*Scrit,dts[j]*2.*np.pi,e0) for i in range(N) for j in range(N)]

pool = InterruptiblePool()
res = pool.map(simulation,parameters)
res = np.nan_to_num(res)
megno = np.clip(res[:,0].reshape((N,N)),1.8,4.)
lyaptime = np.clip(np.absolute(res[:,1].reshape((N,N))),1.,1.e5)/2./np.pi # divide by 2pi to get in units of orbital period

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

f,axarr = plt.subplots(2)
extent=[dts.min(), dts.max(), Ss.min(), Ss.max()]

for ax in axarr:
	ax.set_xlim(extent[0], extent[1])
	ax.set_ylim(extent[2], extent[3])
	ax.set_xlabel(r"$\Delta t / t_{orb}$")
	ax.set_ylabel(r"$S/S_{crit}$")

im1 = axarr[0].imshow(megno, vmin=1.8, vmax=4., aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
cb1 = plt.colorbar(im1, ax=axarr[0])
cb1.solids.set_rasterized(True)
cb1.set_label("MEGNO $\\langle Y \\rangle$")

im2 = axarr[1].imshow(lyaptime, vmin=1., vmax=lyaptime.max(), norm=LogNorm(), aspect='auto', origin="lower", interpolation="nearest", cmap="RdYlGn", extent=extent)
cb2 = plt.colorbar(im2, ax=axarr[1])
cb2.solids.set_rasterized(True)
cb2.set_label("Lyapunov timescale (orbital periods)")

plt.savefig("stark.pdf")

import os
os.system("open stark.pdf")
