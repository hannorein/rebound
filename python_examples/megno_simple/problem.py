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
        sim = rebound.Simulation()
        sim.integrator = integrator
        sim.dt = 0.2
            
        sim.add(m=1.)
        sim.add(m=0.01, a=1,e=0.1)
        sim.add(m=0.01, a=2.)

        sim.move_to_com()
        sim.init_megno()
        particles = sim.particles
        times = np.logspace(2,5,num=1000)
        for t in times:
            sim.integrate(t,0)
            print("%e %e %e %e %e %e %e %e\n" %(sim.t, sim.megno(), particles[0].x, particles[1].x, particles[2].x, particles[3].x, particles[4].x, particles[5].x),file=f)

simulation("whfast")
simulation("ias15")
