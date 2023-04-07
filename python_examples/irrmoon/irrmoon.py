"""
This is an example to investigate the stability of irregular moons. 
This simulation involves a jupiter mass planet, several moons and a star acted as perturber. The moons have different inclinations. I also include some retrograde moons.

The moons undergo Lidov–Kozai oscillation, so some moons will unbound from the planet because of the high eccentricity. The retrograded ones are more bounded.
"""
import rebound
import numpy as np
import matplotlib.pyplot as plt

sim = rebound.Simulation()

#add a jupiter mass planet in the center
jup = rebound.Particle(m=1e-3)
sim.add(jup)

#add some inclined moons
nmoon = 8
moonvar = np.linspace(30,150,nmoon)
moonx = np.linspace(0.2,0.201,nmoon) # slightly different semimajor axis to avoid collision
for i in range(len(moonvar)):
   moon = rebound.Particle()
   incmoon = moonvar[i]
   moon.m = 1e-30 # The moons have small mass to avoid they influence each others
   moon.x = moonx[i]
   moon.vy = np.cos(incmoon/180.*np.pi)*np.sqrt((jup.m+moon.m)/moon.x)
   moon.vz = np.sin(incmoon/180.*np.pi)*np.sqrt((jup.m+moon.m)/moon.x)
   sim.add(moon)

# add a star 
star = rebound.Particle()
star.m = 1.
star.x = 5.
star.vy = np.sqrt((jup.m+star.m)/star.x)
sim.add(star)

sim.move_to_hel()
   
tmax = 1.6e3
sim.integrator = "ias15"
sim.dt = 1e-2

# set up simulation
nout = 1000
xs = np.zeros((nmoon,nout))
eccs = np.zeros((nmoon,nout))
incs = np.zeros((nmoon,nout))
times = np.linspace(0.,tmax,nout)
cons = np.zeros((nmoon,nout))
ps = sim.particles
# record data for each time
for i in range(nout):
    sim.integrate(times[i])
    sim.move_to_hel()
    os = sim.calculate_orbits()
    for j in range(nmoon):
       tar = ps[j+1]
       xs[j][i] = np.sqrt(tar.x**2+tar.y**2+tar.z**2)
       eccs[j][i] = os[j].e
       #incs[j][i] = np.arctan(tar.vz/np.sqrt(tar.vx**2+tar.vy**2))
       #cons[j][i] = np.sqrt(1-tar.e**2)*np.sqrt(tar.vx**2+tar.vy**2)/np.sqrt(tar.vx**2+tar.vy**2+tar.vz**2)
    print(i)
    #print(xs[i])
    #print(eccs[i])
    #print(incs[i])

# plot r from jupiter
for j in range(nmoon):
   plt.plot(times/(2*np.pi/(np.sqrt(10**(-3)/0.2**3))),xs[j],label = 'a = {:.3}'.format(moonvar[j]))
plt.ylabel('r')
plt.xlabel(r't ($P_{moon}$)')
plt.yscale('log')
plt.ylim(bottom = 0.01)
plt.legend()
plt.show()

'''
for j in range(nmoon):
   plt.plot(times/(2*np.pi/(np.sqrt(10**(-3)/0.2**3))),incs[j],label = 'a = {:.3}'.format(moonvar[j]))
plt.ylabel('inc')
plt.xlabel(r't ($P_{moon}$)')
plt.legend()
plt.show()
'''
for j in range(nmoon):
   plt.plot(times/(2*np.pi/(np.sqrt(10**(-3)/0.2**3))),eccs[j],label = 'a = {:.3}'.format(moonvar[j]))
plt.ylabel('ecc')
plt.xlabel(r't ($P_{moon}$)')
plt.yscale('log')
plt.ylim(bottom = 0.01)
plt.legend()
plt.show()

#plt.plot(times,cons[-1])
#plt.show()
