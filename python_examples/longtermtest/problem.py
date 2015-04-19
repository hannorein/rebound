# Import the rebound module
import sys; sys.path.append('../../python_modules')
import rebound
import numpy as np
import time
from interruptible_pool import InterruptiblePool

def simulation(par):
    integrator, run, trial = par
    rebound.reset()
    k = 0.01720209895    
    G = k*k
    rebound.set_G(G)     
    rebound.set_dt(dt)
    if integrator=="mikkola-cor3":
        integrator="mikkola"
        rebound.set_integrator_mikkola_corrector(3)
    elif integrator=="mikkola-cor5":
        integrator="mikkola"
        rebound.set_integrator_mikkola_corrector(5)
    elif integrator=="mikkola-cor7":
        integrator="mikkola"
        rebound.set_integrator_mikkola_corrector(7)
    elif integrator=="mikkola-cor11":
        integrator="mikkola"
        rebound.set_integrator_mikkola_corrector(11)
    else:
        rebound.set_integrator_mikkola_corrector(0)
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
        if integrator=="wh" or integrator=="mercury":
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

    times = np.logspace(np.log10(100.*dt),np.log10(tmax),1000)
    if integrator=="wh" or integrator=="mercury":
        move_to_heliocentric()
    else:
        rebound.move_to_center_of_momentum()
    ei = energy()

    es = []

    runtime = 0.
    for t in times:
        rebound.integrate(t,exactFinishTime=0)
        ef = energy()
        e = np.fabs((ei-ef)/ei)+1.1e-16
        es.append(e)
        runtime += rebound.get_timing()
    
    es = np.array(es)
    print integrator + " done. %.5fs"%(runtime)
    return [times, es]

#3dt = 100.23
dt = 1.3
tmax = 365.*11.8618*1e2
integrators = ["wh","mikkola","ias15","mikkola-cor3","mikkola-cor5","mikkola-cor7","mikkola-cor-11","mercury"]
colors = ["b","r","g","y","m","c","k"]
trials = 4
    
parameters = [(inte,i*trials+j,j) for i,inte in enumerate(integrators) for j in xrange(trials)]
print "Running %d simulations" % (len(parameters))

pool = InterruptiblePool()
res = np.array(pool.map(simulation,parameters)).reshape(len(integrators),trials,2,1000)
print res.shape

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm


f,axarr = plt.subplots(1,1,figsize=(7,5))
extent=[res[:,:,0,:].min()/365./11.8618, res[:,:,0,:].max()/365./11.8618, 1e-16, 1e-5]

axarr.set_xlim(extent[0], extent[1])
axarr.set_ylim(extent[2], extent[3])
axarr.set_xlabel(r"time [Jupiter years]")
axarr.set_ylabel(r"rel energy error")
plt.xscale('log', nonposy='clip')
plt.yscale('log', nonposy='clip')


res_mean = np.mean(res,axis=1)
for i in xrange(len(res)):
    for j in xrange(trials):
        res_trial = res[i,j,:,:]
        im1 = axarr.plot(res_trial[0]/365./11.8618,res_trial[1], color=colors[i],alpha=0.1)
    im1 = axarr.plot(res_mean[i][0]/365./11.8618,res_mean[i][1], label=integrators[i],color=colors[i])

plt.legend(loc='upper left')
plt.savefig("longtermtest.pdf")
import os
os.system("open longtermtest.pdf")
