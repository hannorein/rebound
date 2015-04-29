# Import the rebound module
import sys; sys.path.append('../../python_modules')
import rebound
import numpy as np
import time
from rebound.interruptible_pool import InterruptiblePool

print(rebound.get_build_str())
def simulation(par):
    integrator, run, trial = par
    rebound.reset()
    rebound.set_dt(dt)
    rebound.set_integrator(integrator)
    rebound.set_force_is_velocitydependent(0)

    mfac = 1./1.988544e30
    vfac = 1./0.01720209895    

    rebound.add( m=mfac*1.988544e30, x=3.256101656448802e-03  , y=-1.951205394420489e-04 , z=-1.478264728548705e-04,  vx=vfac*3.039963463108432e-06 ,  vy=vfac*6.030576499910942e-06 ,  vz=vfac*-7.992931269075703e-08, )
    rebound.add( m=mfac*3.302e23, x=-1.927589645545195e-01 , y=2.588788361485397e-01  , z=3.900432597062033e-02 , vx=vfac*-2.811550184725887e-02,  vy=vfac*-1.586532995282261e-02,  vz=vfac*1.282829413699522e-03 , )
    rebound.add( m=mfac*48.685e23, x=-5.976537074581466e-01 , y=3.918678996109574e-01  , z=3.990356741282203e-02 , vx=vfac*-1.113090630745269e-02,  vy=vfac*-1.703310700277280e-02,  vz=vfac*4.089082927733997e-04 ,)
    rebound.add( m=mfac*6.0477246e24, x=-7.986189029000561e-01 , y=-6.086873314992410e-01 , z=-1.250824315650566e-04, vx=vfac*1.012305635253317e-02 ,  vy=vfac*-1.376389620972473e-02,  vz=vfac*3.482505080431706e-07 , )
    rebound.add( m=mfac*6.4185e23, x=7.897942807177173e-01  , y=1.266671734964037e+00  , z=7.092292179885432e-03 , vx=vfac*-1.135279609707971e-02,  vy=vfac*8.579013475676980e-03 ,  vz=vfac*4.582774369441005e-04 , )
    rebound.add( m=mfac*1898.13e24, x=-4.314503046344270e+00 , y=3.168094294126697e+00  , z=8.331048545353310e-02 , vx=vfac*-4.555986691913995e-03,  vy=vfac*-5.727124269621595e-03,  vz=vfac*1.257262404884127e-04 , )
    rebound.add( m=mfac*5.68319e26, x=-4.882304833383455e+00 , y=-8.689263067189865e+00 , z=3.453930436208210e-01 , vx=vfac*4.559352462922572e-03 ,  vy=vfac*-2.748632232963112e-03,  vz=vfac*-1.337915989241807e-04, )
    rebound.add( m=mfac*86.8103e24, x=1.917757033372740e+01  , y=5.671738750949031e+00  , z=-2.273858614425555e-01,  vx=vfac*-1.144087185031310e-03,  vy=vfac*3.588282323722787e-03 ,  vz=vfac*2.829006644043203e-05 , )
    rebound.add( m=mfac*102.41e24, x=2.767031517959636e+01  , y=-1.150331645280942e+01 , z=-4.008018419157927e-01,  vx=vfac*1.183702780101068e-03 ,  vy=vfac*2.917115980784960e-03 ,  vz=vfac*-8.714411604869349e-05, )
    rebound.add( m=mfac*1.4639248e+22, x=7.765250227278298e+00  , y=-3.190996242617413e+01 , z=1.168394015703735e+00 , vx=vfac*3.112825364672655e-03 ,  vy=vfac*1.004673400082409e-04 ,  vz=vfac*-9.111652976208292e-04,)
	
    
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
                E_pot -= particles[i].m*particles[j].m/np.sqrt(r2)
        return E_kin+E_pot

    times = np.logspace(np.log10(100.*dt),np.log10(tmax),Ngrid)
    if integrator=="wh" or integrator=="mercury" or integrator[0:7]=="swifter":
        move_to_heliocentric()
    else:
        rebound.move_to_center_of_momentum()
    ei = energy()

    es = []

    runtime = 0.
    for t in times:
        rebound.integrate(t,exactFinishTime=0,keepSynchronized=0)
        ef = energy()
        e = np.fabs((ei-ef)/ei)+1.1e-16
        es.append(e)
        runtime += rebound.get_timing()
    
    integrator, run, trial = par
    print integrator.ljust(13) + " %9.5fs"%(runtime) + "\t Error: %e"  %( e)
    
    es = np.array(es)
    return [times, es]

Ngrid = 500
#3dt = 100.23
orbit = 11.8618*2.*np.pi
dt = orbit/3000.
tmax = orbit*2e2
integrators = ["mercury","mikkola","wh","mikkola-cor11","swifter-whm"]
#integrators = ["mercury","ias15","wh","mikkola","mikkola-cor3","mikkola-cor5","mikkola-cor7","mikkola-cor11"]
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
trials = 4
    
parameters = [(inte,i*trials+j,j) for i,inte in enumerate(integrators) for j in xrange(trials)]


if len(sys.argv)!=2:
    pool = InterruptiblePool()
    print "Running %d simulations" % (len(parameters))
    res = np.array(pool.map(simulation,parameters)).reshape(len(integrators),trials,2,Ngrid)
    np.save("res.npy",res)
else:
    print "Loading %d simulations" % (len(parameters))
    print sys.argv[1]
    res = np.load(sys.argv[1])

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm


f,axarr = plt.subplots(1,1,figsize=(13,4))
extent=[res[:,:,0,:].min()/orbit, res[:,:,0,:].max()/orbit, 1e-16, 1e-5]

axarr.set_xlim(extent[0], extent[1])
axarr.set_ylim(extent[2], extent[3])
axarr.set_xlabel(r"time [orbits]")
axarr.set_ylabel(r"relative energy error")
plt.xscale('log', nonposy='clip')
plt.yscale('log', nonposy='clip')
plt.grid(True)


res_mean = np.mean(res,axis=1)
for i in xrange(len(res)):
    for j in xrange(trials):
        res_trial = res[i,j,:,:]
        im1 = axarr.plot(res_trial[0]/orbit,res_trial[1], color=colors[integrators[i]],alpha=0.1)
    im1 = axarr.plot(res_mean[i][0]/orbit,res_mean[i][1], label=integrators[i],color=colors[integrators[i]])

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
lgd = plt.legend(loc="upper center",  bbox_to_anchor=(0.5, -0.2),  prop = fontP,ncol=4,frameon=False, numpoints=1, scatterpoints=1 , handletextpad = 0.2, markerscale=2.)
plt.savefig("longtermtest.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
import os
os.system("open longtermtest.pdf")
