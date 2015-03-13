# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool
import time

def energy():
	kin = 0.
	pot = 0.

	N = rebound.get_N()
	N_megno = rebound.get_N_megno()
	particles = rebound.particles_get()
	G = rebound.get_G()

	p = particles[1]
	star = particles[0]
	
	kin = 0.5*(	p.vx**2 + p.vy**2 + p.vz**2)
	pot = -G*star.m/np.sqrt((p.x-star.x)**2 + (p.y-star.y)**2 + (p.z-star.z)**2)
	
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
	rebound.particle_add(primary=sun,m=0.,a=1.,anom=np.pi,e=e0)
	
	rebound.move_to_center_of_momentum()
	rebound.megno_init(1.e-16)
	
	particles = rebound.particles_get()

	def starkforce(): # need to put inside simulation(par) to have access to S and particles
		particles[1].ax += S
	
	rebound.set_additional_forces(starkforce)
	
	#rebound.integrate(50*np.pi)

	E0 = energy()
	
	xs = []
	ys = []	
	E_error = []
	ener = []
	ts = []

	steps = 0
	while rebound.get_t()<1000*np.pi:
		rebound.step()
		steps += 1
		#print(particles[1].x, particles[1].y, particles[0].x, particles[0].y)
		if steps % 1 == 0:
			ts.append(rebound.get_t())
			ener.append(energy())
			rel_E_error = (energy() - E0)/np.fabs(E0) # > 0 if energy increased (came closer to unbound E = 0)
			signed_log_E_error = np.log10(np.fabs(rel_E_error)) if rel_E_error > 0 else -np.log10(np.fabs(rel_E_error)) # hold log, but with sign
			E_error.append(rel_E_error)#signed_log_E_error)
			if np.sqrt(particles[1].x**2 + particles[1].y**2) > 10:
				return [rebound.get_megno(), ts, E_error, ener]
			#xs.append(particles[1].x)
			#ys.append(particles[1].y)
	'''	
	n=4000
	fig, ax = plt.subplots()
	ax.plot(xs[::],ys[::])
	plt.show()
	'''
	
	#return [rebound.get_megno(), 1./(rebound.get_lyapunov()*2.*np.pi)]
	return [rebound.get_megno(), ts, E_error, ener]
#I always set the (osculating) semimajor axis to 1, you can pass different initial e values

e0 = 0.99999999999
Scrit = 0.25 # always true if you use G=M=a=1

N = 100

start_time = time.time()
res = simulation((0.*Scrit,0.01*2*np.pi,e0))
print("Took {0} seconds".format(time.time()-start_time))

ts = [i/2./np.pi for i in res[1]]
E_error = res[2]
energy = res[3]

print("t = ", ts[-1])
print("MEGNO = ", res[0])

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt

f,axarr = plt.subplots(2)

axarr[0].plot(ts, E_error, 'b.')
axarr[0].set_xlabel("Time (orbits)")
axarr[0].set_ylabel("Relative energy error")
'''
axarr[1].plot(ts,energy, 'b.')
axarr[1].set_xlabel("Time (orbits)")
axarr[1].set_ylabel("Energy")
'''
plt.savefig("energyerror.pdf")

import os
os.system("open energyerror.pdf")
