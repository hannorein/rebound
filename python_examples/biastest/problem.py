# Import the rebound module
import sys; sys.path.append('../../python_modules')
import rebound
import numpy as np

def simulation(par):
    integrator, mass = par
    rebound.reset()
    mass = pow(10.,mass)
    k = 0.01720209895    
    G = k*k
    rebound.G = G     
    rebound.dt = 0.
    rebound.integrator = integrator

    rebound.add(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9)   # Sun
    rebound.add(m=mass,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3, vy=+5.51815399480116e-3, vz=-2.66711392865591e-6)   # Jupiter
    rebound.add(m=mass,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3, vy=+3.99723751748116e-3, vz=+1.67206320571441e-5)   # Saturn
    rebound.add(m=mass,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3, vy=+2.06438412905916e-3, vz=-2.17699042180559e-5)   # Uranus
    rebound.add(m=mass,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4, vy=-3.11361111025884e-3, vz=+3.58344705491441e-5)   # Neptune
    N = rebound.N

    def move_to_heliocentric():
        particles = rebound.particles
        
        for i in xrange(1,N):
            particles[i].x -= particles[0].x
            particles[i].y -= particles[0].y
            particles[i].z -= particles[0].z
            particles[i].vx -= particles[0].vx
            particles[i].vy -= particles[0].vy
            particles[i].vz -= particles[0].vz
        particles[0].x  = 0.
        particles[0].y  = 0. 
        particles[0].z  = 0. 
        particles[0].vx = 0. 
        particles[0].vy = 0. 
        particles[0].vz = 0. 


    def energy():
        if integrator=="wh":
            rebound.move_to_com()
        particles = rebound.particles
        E_kin = 0.
        E_pot = 0.
        for i in xrange(N):
            E_kin += 0.5*particles[i].m*(particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz)
            for j in xrange(i+1,N):
                dx = particles[i].x-particles[j].x
                dy = particles[i].y-particles[j].y
                dz = particles[i].z-particles[j].z
                r2 = dx*dx + dy*dy + dz*dz
                E_pot -= G*particles[i].m*particles[j].m/np.sqrt(r2)
        if integrator=="wh":
            move_to_heliocentric()
        return E_kin+E_pot

    rebound.move_to_com()
    ei = energy()

    es = 1e-20

    for s in xrange(1000):
        rebound.step()
        ef = energy()
        e = np.fabs((ei-ef)/ei)
        es = max(es,e)

    return es
N=54
masses = np.linspace(-10.,0.,N)
integrators = ["wh","whfast"]
parameters = [(i,m) for i in integrators for m in masses]

pool = rebound.InterruptiblePool()
res = pool.map(simulation,parameters)
res = np.array(res).reshape(len(integrators),N)


import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm


f,axarr = plt.subplots(1,1,figsize=(7,5))
masses = np.power(10.,masses)
extent=[masses.min(), masses.max(), 9e-21, 1e-1]

axarr.set_xlim(extent[0], extent[1])
axarr.set_ylim(extent[2], extent[3])
axarr.set_xlabel(r"mass ratio")
axarr.set_ylabel(r"rel energy error")
plt.xscale('log', nonposy='clip')
plt.yscale('log', nonposy='clip')


for i in xrange(len(res)):
    print(res[i])
    im1 = axarr.plot(masses,res[i], label=integrators[i])

plt.legend(loc='upper left')
plt.savefig("longtermtest.pdf")
from sys import platform as _platform
if _platform == "darwin":
    import os
    os.system("open longtermtest.pdf")
