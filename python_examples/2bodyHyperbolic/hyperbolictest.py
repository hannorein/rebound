# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool

rebound.reset()
rebound.set_integrator("mikkola")
rebound.set_dt(0.01*2.*np.pi)

rebound.particle_add(m=1.)
rebound.particle_add(m=0., a=1., e=1.1, anom=-10.,MEAN=True)
particles = rebound.particles_get()

for i in xrange(100):
    print particles[1].x, particles[1].y
    rebound.step()
    
