# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
import matplotlib.pyplot as plt
import 
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
	
	#rebound.integrate(1.*2*np.pi)

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
	
	return [rebound.get_megno(), 1./(rebound.get_lyapunov()*2.*np.pi)]

#I always set the (osculating) semimajor axis to 1, you can pass different initial e values

e0 = 0.8 # Rauch uses 0.9 for Fig 4
Scrit = 0.25 # always true if you use G=M=a=1

N = 4
dts = np.linspace(0,2*np.pi,N)
Ss = np.linspace(0,0.5,N)*Scrit
parameters = zip(Ss, dts, [e0]*N)

pool = Interruptiblepool.Pool()
