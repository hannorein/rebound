# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool

def energy():
	kin = 0.
	pot = 0.

	N = rebound.get_N()
	N_megno = rebound.get_N_megno()
	particles = rebound.particles_get()
	G = rebound.get_G()

	for i in range(N-N_megno):
		p = particles[i]
		kin += 0.5*p.m * (p.vx**2 + p.vy**2 + p.vz**2)
		for j in range(i+1,N-N_megno):
			pj = particles[j]
			dx = p.x - pj.x
			dy = p.y - pj.y
			dz = p.z - pj.z
			pot -= G*pj.m*p.m/np.sqrt(dx**2 + dy**2 + dz**2)
	
	return kin + pot

def simulation(par):
	S, dt,e0 = par
	
	rebound.reset()
	rebound.set_integrator("mikkola")
	rebound.set_dt(dt)

	r0 = 1-e0 # assumes the orbit has oculating a = 1 and we start at pericenter
	vy = np.sqrt(2/(1-e0)-1.)
	
	sun = Particle(m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
	rebound.particle_add(sun)
	rebound.particle_add(primary=sun,m=0.,a=1.,anom=0.,e=e0,omega=0.,inc=0.,Omega=0.,MEAN=True)
	
	rebound.move_to_center_of_momentum()
	rebound.megno_init(1.e-16)
	
	particles = rebound.particles_get()

	def starkforce(): # need to put inside simulation(par) to have access to S and particles
		particles[1].ax += S
	
	rebound.set_additional_forces(starkforce)
	
	#rebound.integrate(50*np.pi)

	xs = []
	ys = []	
	ener = []
	steps = 0
	while rebound.get_t()<14*np.pi:
		rebound.step()
		steps += 1
		if steps % 1 == 0:
			ener.append(energy())
			if np.sqrt(particles[1].x**2 + particles[1].y**2) > 10:
				return [rebound.get_megno(), rebound.get_t(), ener]
			#xs.append(particles[1].x)
			#ys.append(particles[1].y)
	'''	
	n=4000
	fig, ax = plt.subplots()
	ax.plot(xs[::],ys[::])
	plt.show()
	'''
	
	#return [rebound.get_megno(), 1./(rebound.get_lyapunov()*2.*np.pi)]
	return [rebound.get_megno(), rebound.get_t(), ener]
#I always set the (osculating) semimajor axis to 1, you can pass different initial e values

e0 = 0.9 # Rauch uses 0.9 for Fig 4
Scrit = 0.25 # always true if you use G=M=a=1

N = 100

res = simulation((0.1*Scrit, 0.18, 0.9))

print(res[1]/2./2/np.pi)
print(res[0])

'''
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
'''
