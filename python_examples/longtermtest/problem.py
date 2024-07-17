# Import the rebound module
import sys
import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
import rebound
import numpy as np
import time
import warnings

def simulation(par):
    integrator, run, trial = par
    sim = rebound.Simulation()
    k = 0.01720209895    
    Gfac = 1./k
    sim.dt = dt
    if integrator == "whfast-nocor":
        integrator = "whfast"
    else:
        sim.ri_whfast.corrector = 11
    sim.integrator = integrator 
    sim.ri_whfast.safe_mode = 0

    massfac = 1.
    sim.add(m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6,      vx=+6.69048890636161e-6*Gfac, vy=-6.33922479583593e-6*Gfac, vz=-3.13202145590767e-9*Gfac)   # Sun
    sim.add(m=massfac/1407.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3*Gfac, vy=+5.51815399480116e-3*Gfac, vz=-2.66711392865591e-6*Gfac)   # Jupiter
    sim.add(m=massfac/3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3*Gfac, vy=+3.99723751748116e-3*Gfac, vz=+1.67206320571441e-5*Gfac)   # Saturn
    sim.add(m=massfac/22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3*Gfac, vy=+2.06438412905916e-3*Gfac, vz=-2.17699042180559e-5*Gfac)   # Uranus
    sim.add(m=massfac/19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4*Gfac, vy=-3.11361111025884e-3*Gfac, vz=+3.58344705491441e-5*Gfac)   # Neptune
    N = sim.N
    particles = sim.particles
    np.random.seed(run)
    for p in particles:
        p.m *= 1.+1e-3*np.random.rand()
        p.x *= 1.+1e-3*np.random.rand()
        p.y *= 1.+1e-3*np.random.rand()
        p.z *= 1.+1e-3*np.random.rand()
        p.vx *= 1.+1e-3*np.random.rand()
        p.vy *= 1.+1e-3*np.random.rand()
        p.vz *= 1.+1e-3*np.random.rand()

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
        if integrator=="mercury" or integrator[0:7]=="swifter":
            mtot = 0.
            for p in particles:
                com_vx += p.vx*p.m 
                com_vy += p.vy*p.m 
                com_vz += p.vz*p.m 
                mtot += p.m
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

    times = np.logspace(np.log10(orbit),np.log10(tmax),Ngrid)
    if integrator=="mercury" or integrator[0:7]=="swifter":
        move_to_heliocentric()
    else:
        sim.move_to_com()
    ei = energy()

    es = []

    runtime = 0.
    start = time.time()
    # Capture warning messages (WHFast timestep too large)
    with warnings.catch_warnings(record=True) as w: 
        warnings.simplefilter("always")
        for t in times:
            sim.integrate(t,exact_finish_time=0)
            ef = energy()
            e = np.fabs((ei-ef)/ei)+1.1e-16
            es.append(e)
    
    integrator, run, trial = par
    print(integrator.ljust(13) + " %9.5fs"%(time.time()-start) + "\t Error: %e"  %( e))
    
    es = np.array(es)
    return [times, es]

Ngrid = 500
#3dt = 100.23
orbit = 11.8618*1.*np.pi
dt = orbit/3000.
tmax = orbit*1e2        # Maximum integration time.
integrators = ["whfast-nocor", "whfast"]
#integrators = ["mercury","swifter-whm","whfast-nocor", "whfast"]
colors = {
    'whfast-nocor': "#FF0000",
    'whfast':       "#00AA00",
    'mercury':      "#6E6E6E",
    'swifter-whm':  "#444444",
    'swifter-helio':"#AABBBB",
    'swifter-tu4':  "#FFAAAA",
    'ias15':        "g",
}
trials = 4
    
parameters = [(inte,i*trials+j,j) for i,inte in enumerate(integrators) for j in xrange(trials)]
if len(sys.argv)!=2:
    try:
        from multiprocess import Pool
    except:
        raise RuntimeError("Please install the multiprocess module with `pip install multiprocess`.")
    with Pool() as pool:
        print("Running %d simulations on %d threads..." % (len(parameters), pool._processes))
        res = np.array(pool.map(simulation,parameters)).reshape(len(integrators),trials,2,Ngrid)
    np.save("res.npy",res)
else:
    print("Loading %d simulations" % (len(parameters)))
    print(sys.argv[1])
    res = np.load(sys.argv[1])



f,axarr = plt.subplots(1,1,figsize=(13,4))
extent=[res[:,:,0,:].min()/orbit, res[:,:,0,:].max()/orbit, 1e-16, 1e-5]

axarr.set_xlim(extent[0], extent[1])
axarr.set_ylim(extent[2], extent[3])
axarr.set_xlabel(r"time [orbits]")
axarr.set_ylabel(r"relative energy error")
plt.xscale('log')
plt.yscale('log')
plt.grid(True)


res_mean = np.mean(res,axis=1)
for i in xrange(len(res)):
    for j in xrange(trials):
        res_trial = res[i,j,:,:]
        im1 = axarr.plot(res_trial[0]/orbit,res_trial[1], color=colors[integrators[i]],alpha=0.2)
    im1 = axarr.plot(res_mean[i][0]/orbit,res_mean[i][1], label=integrators[i].upper(),color=colors[integrators[i]], linewidth=2.0)

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
lgd = plt.legend(loc="upper center",  bbox_to_anchor=(0.5, -0.2),  prop = fontP,ncol=3,frameon=False, numpoints=1, scatterpoints=1 , handletextpad = 0.2, markerscale=2.)
plt.savefig("longtermtest.pdf", bbox_extra_artists=(lgd,), bbox_inches='tight')
from sys import platform as _platform
if _platform == "darwin":
    import os
    os.system("open longtermtest.pdf")
