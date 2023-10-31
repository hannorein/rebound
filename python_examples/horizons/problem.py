import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
import rebound
import socket
import sys
import os.path
import os
filename = "cache.bin"

solar_system_objects = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "C/2014 Q2"]

if os.path.isfile(filename):
    # Try to load simulation from file
    sim = rebound.Simulation(filename)
else: 
    sim = rebound.Simulation()
    # Get data from NASA Horizons
    try:
        sim.add(solar_system_objects)
    except socket.error:
        print("A socket error occured. Maybe Horizons is down?")
        sys.exit(0) # we ignore the error and exit

    sim.move_to_com()
    # Configure simulation
    sim.integrator = "whfast"
    sim.set_dt = 0.01
    # Let's save it for next time
    # Note: sim.save_to_file() only saves the particle data, not the integrator settings, etc.
    sim.save_to_file(filename)  

sim.status()

import numpy as np
Nout = 1000
times = np.linspace(0,16.*np.pi,Nout) # 8 years
x = np.zeros((sim.N,Nout))
y = np.zeros((sim.N,Nout))

ps = sim.particles
for ti,t in enumerate(times):
    sim.integrate(t)
    for i, p in enumerate(ps):
        x[i][ti] = p.x
        y[i][ti] = p.y


fig = plt.figure(figsize=(11,5))

def plot(zoom):
    ax.set_xlim([-zoom,zoom])
    ax.set_ylim([-zoom,zoom])
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    for i in xrange(0,sim.N):
        plt.plot(x[i],y[i])
        if x[i][-1]*x[i][-1]+y[i][-1]*y[i][-1]>0.01*zoom*zoom or i==0:
            ax.annotate(solar_system_objects[i], xy=(x[i][-1], y[i][-1]),horizontalalignment="center")

ax = plt.subplot(121)
plot(zoom=24.)
ax = plt.subplot(122)
plot(zoom=1.2)

plt.savefig("orbits.pdf")
