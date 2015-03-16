# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool


for i in np.linspace(-2.*np.pi,2.*np.pi,1000):
    rebound.reset()
    rebound.set_integrator("mikkola")
    rebound.set_dt(0.01*2.*np.pi)

    try:
        rebound.particle_add(m=1.)
        rebound.particle_add(m=0., a=1., e=1.01, anom=i)
        particles = rebound.particles_get()
        print particles[1].x, particles[1].y, i
        #rebound.step()
        #print particles[1].x, particles[1].y, i
    except:
        pass
    
