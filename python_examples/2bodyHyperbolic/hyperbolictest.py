# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle
import numpy as np
from interruptible_pool import InterruptiblePool


for i in np.linspace(-1.67,1.67,100):
    rebound.reset()
    rebound.set_integrator("mikkola")
    rebound.set_dt(0.0002*2.*np.pi)

    try:
        rebound.particle_add(m=1.)
        rebound.particle_add(m=0., a=1., e=11, anom=i)
        particles = rebound.particles_get()
        print particles[1].x, particles[1].y
        rebound.step()
        print particles[1].x, particles[1].y
    except:
        pass
    
