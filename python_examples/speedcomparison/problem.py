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
    rebound.G = G     
    rebound.dt = dt
    if integrator == "whfast-nocor":
        integrator = "whfast"
    else:
        rebound.integrator_whfast_corrector = 11
    rebound.integrator = integrator
    rebound.force_is_velocitydependent = 0

    rebound.add(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9)   # Sun
    rebound.add(m=1./1407.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3, vy=+5.51815399480116e-3, vz=-2.66711392865591e-6)   # Jupiter
    rebound.add(m=1./3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3, vy=+3.99723751748116e-3, vz=+1.67206320571441e-5)   # Saturn
    rebound.add(m=1./22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3, vy=+2.06438412905916e-3, vz=-2.17699042180559e-5)   # Uranus
    rebound.add(m=1./19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4, vy=-3.11361111025884e-3, vz=+3.58344705491441e-5)   # Neptune
    N = rebound.N
    particles = rebound.particles
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
        particles[0].x  = 0.
        particles[0].y  = 0. 
        particles[0].z  = 0. 
        particles[0].vx = 0. 
        particles[0].vy = 0. 
        particles[0].vz = 0. 


    def energy():
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
        rebound.move_to_com()
    ei = energy()


    runtime = 0.
    rebound.integrate(tmax,exact_finish_time=0)
    ef = energy()
    e = np.fabs((ei-ef)/ei)+1.1e-16
    runtime += rebound.timing
    
    integrator, dt, run = par
    print integrator.ljust(13) + " %9.5fs"%(runtime) + "\t Error: %e"  %( e)
    return [runtime, e]

orbit = 365.*11.8618
dts = np.logspace(-2.5,2,100)
#dts = np.logspace(-3,2,155)
tmax = orbit*1e3
integrators = ["wh","swifter-whm","swifter-tu4","swifter-helio","mercury","whfast-nocor","whfast"]
colors = {
    'whfast-nocor': "#FF0000",
    'whfast':       "#00AA00",
    'mercury':      "#6E6E6E",
    'wh':           "b",
    'swifter-whm':  "#444444",
    'swifter-helio':"#AABBBB",
    'swifter-tu4':  "#FFAAAA",
    'ias15':        "g",
    }
parameters = [(inte,dt,i*len(dts)+j) for i,inte in enumerate(integrators) for j, dt in enumerate(dts)]
   
if len(sys.argv)!=2:
    pool = InterruptiblePool(8)
    print "Running %d simulations" % (len(parameters))
    res = np.array(pool.map(simulation,parameters)).reshape(len(integrators),len(dts),2)
    np.save("res.npy",res)
else:
    print "Loading %d simulations" % (len(parameters))
    print sys.argv[1]
    res = np.load(sys.argv[1])

print res.shape

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')

fig = plt.figure(figsize=(6.5,4))

res_mean = np.mean(res,axis=1)
extent=[1e-16, 1e-5]

ax = plt.subplot(1,1,1)
ax.set_xlim(res[:,:,0].min(), res[:,:,0].max())
ax.set_ylim(extent[0], extent[1])
ax.set_ylabel(r"relative energy error")
ax.set_xlabel(r"runtime [s]")
plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')
plt.grid(True)
for i in xrange(len(integrators)):
    res_i = res[i,:,:]
    im1 = ax.scatter(res_i[:,0], res_i[:,1], label=integrators[i].upper(),color=colors[integrators[i]],s=10)
    #im1 = axarr.scatter(dts, res_i[:,1], label=integrators[i],color=colors[i])


lgd = plt.legend(loc="upper center",  bbox_to_anchor=(0.5, -0.15),  prop = fontP,ncol=3,frameon=False, numpoints=1, scatterpoints=1 , handletextpad = -0.5, markerscale=2.)
plt.savefig("speed.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close(0)

fig = plt.figure(figsize=(6.5,4))
ax = plt.subplot(1,1,1)
ax.set_xlim(orbit/dts.max(), orbit/dts.min())
ax.set_ylim(extent[0], extent[1])
ax.set_ylabel(r"relative energy error")
ax.set_xlabel(r"steps per orbit")
plt.yscale('log', nonposy='clip')
plt.xscale('log', nonposy='clip')
plt.grid(True)
for i in xrange(len(integrators)):
    res_i = res[i,:,:]
    im1 = ax.scatter(orbit/dts, res_i[:,1], label=integrators[i].upper(),color=colors[integrators[i]],s=10)

lgd = plt.legend(loc="upper center",  bbox_to_anchor=(0.5, -0.15),  prop = fontP,ncol=3,frameon=False, numpoints=1, scatterpoints=1 , handletextpad = -0.5, markerscale=2.)

plt.savefig("shorttermenergy.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
