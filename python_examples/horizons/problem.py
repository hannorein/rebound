import rebound
import os.path

filename = "cache.bin"

solar_system_objects = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "C/2014 Q2"]
if os.path.isfile(filename):
    rebound.load(filename)
else: 
    # Get data from NASA Horizons
    rebound.add(solar_system_objects)
    rebound.move_to_com()
    # Let's save it for next time
    # Note: rebound.save() only saves the particle data, not the integrator settings, etc.
    rebound.save("cache.bin")  

rebound.integrator = "whfast"
rebound.set_dt = 0.01
rebound.status()

import numpy as np
Nout = 1000
times = np.linspace(0,16.*np.pi,Nout) # 8 years
x = np.zeros((rebound.N,Nout))
y = np.zeros((rebound.N,Nout))

ps = rebound.particles
for ti,t in enumerate(times):
    rebound.integrate(t)
    for i, p in enumerate(ps):
        x[i][ti] = p.x
        y[i][ti] = p.y

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(11,5))

def plot(zoom):
    ax.set_xlim([-zoom,zoom])
    ax.set_ylim([-zoom,zoom])
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    for i in xrange(0,rebound.N):
        plt.plot(x[i],y[i])
        if x[i][-1]*x[i][-1]+y[i][-1]*y[i][-1]>0.01*zoom*zoom or i==0:
            ax.annotate(solar_system_objects[i], xy=(x[i][-1], y[i][-1]),horizontalalignment="center")

ax = plt.subplot(121)
plot(zoom=24.)
ax = plt.subplot(122)
plot(zoom=1.2)

plt.savefig("orbits.pdf")
