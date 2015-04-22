# Import the rebound module
import numpy as np
import time
import sys
from rebound.interruptible_pool import InterruptiblePool

def simulation(par):
    import rebound
    integrator, dt, run = par
    rebound.reset()
    k = 0.01720209895    
    G = k*k
    rebound.set_G(G)     
    rebound.set_dt(dt)
    rebound.set_integrator(integrator)
    rebound.set_force_is_velocitydependent(0)

    rebound.add_particle(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9)   # Sun
    rebound.add_particle(m=1./1407.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3, vy=+5.51815399480116e-3, vz=-2.66711392865591e-6)   # Jupiter
    rebound.add_particle(m=1./3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3, vy=+3.99723751748116e-3, vz=+1.67206320571441e-5)   # Saturn
    rebound.add_particle(m=1./22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3, vy=+2.06438412905916e-3, vz=-2.17699042180559e-5)   # Uranus
    rebound.add_particle(m=1./19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4, vy=-3.11361111025884e-3, vz=+3.58344705491441e-5)   # Neptune
#    rebound.add_particle(m=0,             x=-2.13858977531573e+1, y=+3.20719104739886e+1, z=+2.49245689556096e+0, vx=-1.76936577252484e-3, vy=-2.06720938381724e-3, vz=+6.58091931493844e-4)   # Pluto
    N = rebound.get_N()
    particles = rebound.get_particles()
    np.random.seed(run)
    for i in xrange(N):
        particles[i].m *= 1.+1e-3*np.random.rand()
        particles[i].x *= 1.+1e-3*np.random.rand()
        particles[i].y *= 1.+1e-3*np.random.rand()
        particles[i].z *= 1.+1e-3*np.random.rand()
        particles[i].vx *= 1.+1e-3*np.random.rand()
        particles[i].vy *= 1.+1e-3*np.random.rand()
        particles[i].vz *= 1.+1e-3*np.random.rand()

    def move_to_heliocentric():
        particles = rebound.get_particles()
        
        particles[0].x  = 0.
        particles[0].y  = 0. 
        particles[0].z  = 0. 
        particles[0].vx = 0. 
        particles[0].vy = 0. 
        particles[0].vz = 0. 


    def energy():
        particles = rebound.get_particles()
        com_vx = 0.
        com_vy = 0.
        com_vz = 0.
        if integrator=="wh" or integrator=="mercury" or integrator[0:7]=="swifter":
            mtot = 0.
            for i in xrange(0,N):
                com_vx += particles[i].vx*particles[i].m 
                com_vy += particles[i].vy*particles[i].m 
                com_vz += particles[i].vz*particles[i].m 
                mtot += particles[i].m
            com_vx /= mtot
            com_vy /= mtot
            com_vz /= mtot
        E_kin = 0.
        E_pot = 0.
        for i in xrange(N):
            dvx = particles[i].vx - com_vx
            dvy = particles[i].vy - com_vy
            dvz = particles[i].vz - com_vz
            E_kin += 0.5*particles[i].m*(dvx*dvx + dvy*dvy + dvz*dvz)
            for j in xrange(i+1,N):
                dx = particles[i].x-particles[j].x
                dy = particles[i].y-particles[j].y
                dz = particles[i].z-particles[j].z
                r2 = dx*dx + dy*dy + dz*dz
                E_pot -= G*particles[i].m*particles[j].m/np.sqrt(r2)
        return E_kin+E_pot

    if integrator=="wh" or integrator=="mercury" or integrator[0:7]=="swifter":
        move_to_heliocentric()
    else:
        rebound.move_to_center_of_momentum()
    ei = energy()


    runtime = 0.
    rebound.integrate(tmax,exactFinishTime=0)
    ef = energy()
    e = np.fabs((ei-ef)/ei)+1.1e-16
    runtime += rebound.get_timing()
    
    integrator, dt, run = par
    print integrator.ljust(13) + " %9.5fs"%(runtime) + "\t Error: %e"  %( e)
    return [runtime, e]

orbit = 365.*11.8618
dts = np.logspace(-2.5,2,128)
#dts = np.logspace(-3,2,155)
tmax = orbit*1e2
integrators = ["wh","mikkola","swifter-whm","swifter-tu4","swifter-helio","mikkola-cor3","mikkola-cor5","mikkola-cor7","mikkola-cor11","mercury"]
colors = {
    'mikkola':      "#FF0000",
    'mikkola-cor3': "#FF7700",
    'mikkola-cor5': "#FF9D00",
    'mikkola-cor7': "#FFC400",
    'mikkola-cor11':"#FFDD00",
    'mikkola-jac':  "#D4FF00",
    'mercury':      "#6E6E6E",
    'wh':           "b",
    'swifter-whm':  "#444444",
    'swifter-helio':"#AABBBB",
    'swifter-tu4':  "#FFAAAA",
    'ias15':        "g",
    }
parameters = [(inte,dt,i*len(dts)+j) for i,inte in enumerate(integrators) for j, dt in enumerate(dts)]
   
if len(sys.argv)!=2:
    pool = InterruptiblePool(12)
    print "Running %d simulations" % (len(parameters))
    res = np.array(pool.map(simulation,parameters)).reshape(len(integrators),len(dts),2)
else:
    print "Loading %d simulations" % (len(parameters))
    print sys.argv[1]
    res = np.load(sys.argv[1])

print res.shape

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm

np.save("res.npy",res)
fig = plt.figure(figsize=(13,4))

res_mean = np.mean(res,axis=1)
extent=[1e-16, 1e-5]

ax = plt.subplot(1,2,1)
ax.set_xlim(res[:,:,0].min(), res[:,:,0].max())
ax.set_ylim(extent[0], extent[1])
ax.set_ylabel(r"relative energy error")
ax.set_xlabel(r"runtime [s]")
plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')
plt.grid(True)
for i in xrange(len(integrators)):
    res_i = res[i,:,:]
    im1 = ax.scatter(res_i[:,0], res_i[:,1], label=integrators[i],color=colors[integrators[i]],s=10)
    #im1 = axarr.scatter(dts, res_i[:,1], label=integrators[i],color=colors[i])

ax = plt.subplot(1,2,2)
ax.set_xlim(orbit/dts.max(), orbit/dts.min())
ax.set_ylim(extent[0], extent[1])
ax.set_ylabel(r"relative energy error")
ax.set_xlabel(r"steps per orbit")
plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')
plt.grid(True)
for i in xrange(len(integrators)):
    res_i = res[i,:,:]
    im1 = ax.scatter(orbit/dts, res_i[:,1], label=integrators[i],color=colors[integrators[i]],s=10)

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
lgd = plt.legend(loc="upper center",  bbox_to_anchor=(-0.1, -0.2),  prop = fontP,ncol=5,frameon=False, numpoints=1, scatterpoints=1 , handletextpad = -0.5, markerscale=2.)
plt.savefig("speed.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
import os
os.system("open speed.pdf")
