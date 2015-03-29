# Import the rebound module
import sys; sys.path.append('../')
import rebound
import numpy as np
from interruptible_pool import InterruptiblePool

def simulation(par):
    integrator = par
    rebound.reset()
    k = 0.01720209895    
    G = k*k
    rebound.set_G(G)     
    rebound.set_dt(dt)
    rebound.set_integrator(integrator)

    rebound.add_particle(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9)   # Sun
    rebound.add_particle(m=1./1047.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3, vy=+5.51815399480116e-3, vz=-2.66711392865591e-6)   # Jupiter
#    rebound.add_particle(m=1./3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3, vy=+3.99723751748116e-3, vz=+1.67206320571441e-5)   # Saturn
#    rebound.add_particle(m=1./22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3, vy=+2.06438412905916e-3, vz=-2.17699042180559e-5)   # Uranus
#    rebound.add_particle(m=1./19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4, vy=-3.11361111025884e-3, vz=+3.58344705491441e-5)   # Neptune
#    rebound.add_particle(m=0,             x=-2.13858977531573e+1, y=+3.20719104739886e+1, z=+2.49245689556096e+0, vx=-1.76936577252484e-3, vy=-2.06720938381724e-3, vz=+6.58091931493844e-4)   # Pluto
    N = rebound.get_N()

    def move_to_heliocentric():
        particles = rebound.get_particles()
        
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
            rebound.move_to_center_of_momentum()
        particles = rebound.get_particles()
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

    times = np.logspace(np.log10(1.*dt),np.log10(tmax),1000)
    rebound.move_to_center_of_momentum()
    ei = energy()

    es = []

    for time in times:
        rebound.integrate(time,exactFinishTime=0)
        ef = energy()
        e = np.fabs((ei-ef)/ei)+1.1e-16
        es.append(e)

    es = np.array(es)
    print integrator + " done."
    return [times, es]

dt = .001
tmax = 365.*1e2
integrators = ["wh","mikkola","ias15"]
    
parameters = [(i) for i in integrators]

pool = InterruptiblePool()
res = pool.map(simulation,parameters)

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm


f,axarr = plt.subplots(1,1,figsize=(7,5))
extent=[res[0][0].min()/365., res[0][0].max()/365., 1e-16, 1e-5]

axarr.set_xlim(extent[0], extent[1])
axarr.set_ylim(extent[2], extent[3])
axarr.set_xlabel(r"time [days]")
axarr.set_ylabel(r"rel energy error")
plt.xscale('log', nonposy='clip')
plt.yscale('log', nonposy='clip')


for i in xrange(len(res)):
    im1 = axarr.plot(res[i][0]/365.,res[i][1], label=integrators[i])

plt.legend(loc='upper left')
plt.savefig("longtermtest.pdf")
import os
os.system("open longtermtest.pdf")
