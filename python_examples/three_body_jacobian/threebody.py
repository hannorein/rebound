__author__ = 'Alysa Obertas'
__email__ = 'obertas@astro.utoronto.ca'

# threebody.py
# circular restricted three body problem
# orbital dynamics of a test particle for a planet with a=1
# the value of mu and the initial conditions of the test
# particle can be changed to see how the trajectory is altered
# Written by Alysa Obertas December 2015

import rebound
from matplotlib import pyplot as plt
import numpy as np

outfile = "three.npz"

###########################################################
## system initialisation

mu = 0.1 # deviation of star mass from 1

mstar = 1-mu # mass of star
mplanet = mu # mass of planet

xtest = 0 # initial x position of test particle
ytest = -1.4
vxtest = 0.99 # initial velocity of test particle in x direction
vytest = 0

# xtest = -1.01 # initial x position of test particle
# ytest = 0
# vxtest = 0 # initial velocity of test particle in x direction
# vytest = -0.96

# xtest = 1.2 # initial x position of test particle
# ytest = 0
# vxtest = 0 # initial velocity of test particle in x direction
# vytest = 1.1

# xtest = 1.3 # initial x position of test particle
# ytest = 0
# vxtest = 0 # initial velocity of test particle in x direction
# vytest = 1.2

# xtest = 0.5 # initial x position of test particle
# ytest = 0
# vxtest = -0.5 # initial velocity of test particle in x direction
# vytest = 0.5

###########################################################
## simulation initialisation

Noutputs= 1e4
Nbody = 2
year = 2.*np.pi # One year in units where G=1
tf = 5 # end time in years
times = np.linspace(0.,tf*year, Noutputs)

x = np.zeros((Nbody,Noutputs))
y = np.zeros((Nbody,Noutputs))

a = np.zeros((Nbody,Noutputs)) # semimajor axes
e = np.zeros((Nbody,Noutputs)) # eccentricities

###########################################################
## start simulation

sim = rebound.Simulation()
sim.add(m=mstar)
sim.add(m=mplanet, a=1)
sim.add(x=xtest, y=ytest, vx=vxtest, vy=vytest)
sim.status()

sim.dt = 0.000001
sim.move_to_com()        
ps = sim.particles       

for i,time in enumerate(times):
    sim.integrate(time)

    a[0][i] = ps[1].a
    e[0][i] = ps[1].e
    a[1][i] = ps[2].a
    e[1][i] = ps[2].e
    
    x[0][i] = ps[1].x
    y[0][i] = ps[1].y
    x[1][i] = ps[2].x
    y[1][i] = ps[2].y

sim.status()

###########################################################
## save time, a, e to file
np.savez(outfile,t=times/year,a=a,e=e,x=x,y=y,vxtest=vxtest,vytest=vytest)




