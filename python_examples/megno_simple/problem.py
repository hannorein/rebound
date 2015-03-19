#!/usr/bin/python
# This example integrates Jupiter and Saturn in the Solar system for a variety of initial conditions.
# Alongside the normal equations of motions, IAS15 is used to integrate the variational equations.
# These can be used to measure the Mean Exponential Growth of Nearby Orbits (MEGNO), a chaos indicator.
# This example script runs 12^2 simulations and plots the MEGNO value. Values close to <Y>=2 correspond 
# to regular quasi-periodic orbits. Higher values of <Y> correspond to chaotic orbits.
 
# Import the rebound module
import sys; sys.path.append('../')
import rebound
# Import other modules
import numpy as np
import multiprocessing

rebound.reset()
rebound.set_integrator("mikkola")
rebound.set_dt(0.1)
    
rebound.particle_add(m=1.)
rebound.particle_add(m=0.0000, a=5.204, anom=0.600, omega=0.257, e=0.048)

rebound.move_to_center_of_momentum()
rebound.megno_init(1e-16)
rebound.integrate(1e4*2.*np.pi)

print(rebound.get_megno())
