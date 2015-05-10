import rebound
import os.path

filename = "cache.bin"
if os.path.isfile(filename):
    rebound.load(filename)
else: 
    rebound.add("Sun")
    rebound.add("Mercury")
    rebound.add("Earth") 
    rebound.add("C/2014 Q2") # Comets!
    rebound.move_to_com()
    rebound.save("cache.bin")

rebound.set_integrator("whfast")
rebound.set_dt(0.01)
rebound.status()

import numpy as np
Nout = 100
times = np.linspace(0,20.*np.pi,Nout) # 1 year
N = rebound.get_N()
x = np.zeros((N,Nout))
y = np.zeros((N,Nout))

p = rebound.get_particles()
for ti,t in enumerate(times):
    rebound.integrate(t)
    for i in xrange(0,N):
        x[i][ti] = p[i].x
        y[i][ti] = p[i].y

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt

for i in xrange(0,N):
    plt.scatter(x[i],y[i])


plt.savefig("orbits.pdf", bbox_inches='tight')
