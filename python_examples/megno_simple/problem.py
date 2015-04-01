#!/usr/bin/python
# This example integrates Jupiter and Saturn in the Solar system for a variety of initial conditions.
# Alongside the normal equations of motions, IAS15 is used to integrate the variational equations.
# These can be used to measure the Mean Exponential Growth of Nearby Orbits (MEGNO), a chaos indicator.
# This example script runs 12^2 simulations and plots the MEGNO value. Values close to <Y>=2 correspond 
# to regular quasi-periodic orbits. Higher values of <Y> correspond to chaotic orbits.
 
# Import the rebound module
import sys; sys.path.append('../../python_modules')
import rebound
from interruptible_pool import InterruptiblePool
# Import other modules
import numpy as np

rebound.reset()
#rebound.set_integrator("ias15")
rebound.set_integrator("mikkola")
rebound.set_dt(0.001)
    
sun = rebound.Particle(m=1)
rebound.add_particle(sun)
planet = rebound.Particle(m=0.01,a=1,e=0.1)
rebound.add_particle(planet)

rebound.add_particle(a=2.)

rebound.move_to_center_of_momentum()
rebound.init_megno(1e-10)
particles = rebound.get_particles()
for x in xrange(100000):
    rebound.step()
    print rebound.get_t(), particles[4].x 

