import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
import rebound
import os.path
import os
import numpy as np

filename = "cache.bin"

# Solar system objects: Sun and the eight planets + Pluto
solar_system_objects = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune","Pluto"]
# Spacecrafts IDs
# You can add other IDs as long as you also add the spacecraft's respective name in the spacecrafts_names list below
spacecrafts_ids = ["-98","-32","-48","-82","-23","-24"]
# Spacecrafts names
spacecrafts_names = ["New Horizon","Voyager 2","HST","Cassini","Pioneer 10","Pioneer 11"]

# creating arrays for all the objects and their corresponding names (the spacecrafts can be found in the HORISONS database based on their IDs rather than their names)
objects = solar_system_objects + spacecrafts_ids
objects_names = solar_system_objects + spacecrafts_names

# We determine the edge of the Solar System as the heliopause
# According to NASA, Voyager 1 reached the heliopause on August 25, 2012
# http://www.nasa.gov/mission_pages/voyager/voyager20130912.html
# This is what we will use to find the edge of the solar system

simv1 = rebound.Simulation()
date = "2012-08-25 10:05"
simv1.add("-31", date=date) # -31 is Voyager 1's ID in HORIZONS
pv1 = simv1.particles
for i,p in enumerate(pv1):
	dx = p.x
	dy = p.y
	dz = p.z

# rv1 is the radial distance in AU that Voyager 1 had reached at the heliopause. We will use it as reference point of where the solar system ends
rv1 = np.sqrt(dx*dx + dy*dy + dz*dz)

if os.path.isfile(filename):
    # Try to load simulation from file
    sim = rebound.Simulation.from_file(filename)
else: 
    sim = rebound.Simulation()
    # Get data from NASA Horizons
    sim.add(objects, date=date)
    sim.move_to_com()
    # Configure simulation
    sim.integrator = "whfast"
    sim.set_dt = 0.01
    # Let's save it for next time
    # Note: sim.save() only saves the particle data, not the integrator settings, etc.
    sim.save(filename)  

sim.status()

#import numpy as np
Nout = 10000
times = np.linspace(0,100.*np.pi,Nout) # 50 years, it should be enough for the latest spacecraft to leave the Solar System, since it took about 35 years for Voyager 1 to do so

x = np.zeros((sim.N,Nout))
y = np.zeros((sim.N,Nout))
z = np.zeros((sim.N,Nout))
r = np.zeros((sim.N,Nout))

ps = sim.particles
for ti,t in enumerate(times):
    sim.integrate(t)
    for i, p in enumerate(ps):
        x[i][ti] = p.x
	#print ti,i,p.x
        y[i][ti] = p.y
        z[i][ti] = p.z
	r[i][ti] = np.sqrt(p.x*p.x + p.y*p.y + p.z*p.z)
	#print r[i][ti]

# Verufy when spacecraft has exited the Solar System
n = len(solar_system_objects) # number of solar system objects
for si,s in enumerate(spacecrafts_names):
	r_spacecraft = r[si+n][:]
	print "Finding info on " + s
	for ti,t in enumerate(times):
		if r_spacecraft[ti] >= rv1:
			print s + " will leave the solar system in " + str(t/(2.0*np.pi)) + " years."
			break
	if r_spacecraft[-1] < rv1:
		print s + " will not leave the Solar Sytem in the next 50 years."

fig = plt.figure(figsize=(11,11))

def plot(zoom): 
    ax.set_xlim([-zoom,zoom])
    ax.set_ylim([-zoom,zoom])
    ax.set_xlabel("x [AU]")
    ax.set_ylabel("y [AU]")
    for i in xrange(0,sim.N):
        plt.plot(x[i],y[i])
        if x[i][-1]*x[i][-1]+y[i][-1]*y[i][-1]>0.01*zoom*zoom or i==0:
		if objects_names[i] == "Cassini": # align the name of Cassini so it does not overlap Saturn's name
			ax.annotate(objects_names[i], xy=(x[i][-1000], y[i][-1000]),horizontalalignment="left")
		else:
           		ax.annotate(objects_names[i], xy=(x[i][-1], y[i][-1]),horizontalalignment="center")

ax = plt.subplot(221)
plot(zoom=1.7) 
ax = plt.subplot(222)
plot(zoom=50.0)
ax = plt.subplot(223)
plot(zoom=250.0)

plt.savefig("orbits.pdf")
