#!/usr/bin/python
# This example script runs 2 simulations and plots the MEGNO value. Values close to <Y>=2 correspond 
# to regular quasi-periodic orbits. Higher values of <Y> correspond to chaotic orbits.
from __future__ import print_function
# Import the rebound module
import rebound
# Import other modules
import numpy as np

def simulation(integrator):
    print("Running "+integrator)
    with open(integrator+".txt","w") as f:
        rebound.reset()
        rebound.integrator = integrator
        rebound.dt = 0.2
            
        rebound.add(m=1.)
        rebound.add(m=0.01, a=1,e=0.1)
        rebound.add(m=0.01, a=2.)

        rebound.move_to_com()
        rebound.init_megno(1e-10)
        particles = rebound.particles
        times = np.logspace(2,5,num=1000)
        for t in times:
            rebound.integrate(t,0)
            print("%e %e %e %e %e %e %e %e\n" %(rebound.t, rebound.calculate_megno(), particles[0].x, particles[1].x, particles[2].x, particles[3].x, particles[4].x, particles[5].x),file=f)

simulation("ias15")
simulation("whfast")
