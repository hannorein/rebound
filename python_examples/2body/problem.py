# Import the rebound module
import rebound
import numpy as np

torb = 2.*np.pi
tmax = 100.34476128*torb

def simulation(par):
    anom, dt, e, integrator = par

    e = 1.-pow(10.,e)
    dt = pow(10.,dt)*torb

    rebound.reset()
    rebound.integrator = integrator
    rebound.integrator_whfast_safe_mode = 0
    rebound.force_is_velocitydependent = 0
    rebound.dt = dt

    rebound.add(m=1.)
    rebound.add(m=0., x=(1.-e), vy=np.sqrt((1.+e)/(1.-e)))
    particles = rebound.particles
    
    Ei = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)
    
    rebound.integrate(tmax,exact_finish_time=0)
    
    Ef = -1./np.sqrt(particles[1].x*particles[1].x+particles[1].y*particles[1].y+particles[1].z*particles[1].z) + 0.5 * (particles[1].vx*particles[1].vx+particles[1].vy*particles[1].vy+particles[1].vz*particles[1].vz)

    return [float(rebound.iter)/rebound.t*dt, np.fabs((Ef-Ei)/Ei)+1e-16, rebound.timing/rebound.t*dt*1e6/2., (Ef-Ei)/Ei]


N = 25
dts = np.linspace(-3.1,-0.1,N)
e0s = np.linspace(0,-8,N)
integrators= ["wh","whfast"]

niter = []
energyerror = []
energyerror_sign = []
timing = []

for integrator in integrators:
    print("Running "+ integrator)
    parameters = [(1.732, dts[i], e0s[j], integrator) for j in range(N) for i in range(N)]

    pool = rebound.InterruptiblePool(2)
    res = pool.map(simulation,parameters)
    res = np.nan_to_num(res)
    niter.append(res[:,0].reshape((N,N)))
    energyerror.append(res[:,1].reshape((N,N)))
    timing.append(res[:,2].reshape((N,N)))
    energyerror_sign.append(res[:,3].reshape((N,N)))

import matplotlib; matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm

f,axarr = plt.subplots(3,2,figsize=(13,14),sharex='col', sharey='row')
extent=[dts.min(), dts.max(), e0s.max(), e0s.min()]
plt.subplots_adjust(wspace = 0.05)
plt.subplots_adjust(hspace = 0.05)

for ay in axarr:
    for ax in ay:
        # ax = axarr
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        if ax.is_last_row():
            ax.set_xlabel(r"timestep, log10$(dt/t_{orb})$")
        if ax.is_first_col():
            ax.set_ylabel(r"eccentricity, log10$(1-e)$")


for i, integrator in enumerate(integrators):
    if axarr[0,i].is_first_row():
        axarr[0,i].set_title(integrator.upper(),fontsize = 16)
    im1 = axarr[0,i].imshow(energyerror[i], norm=LogNorm(), vmax=1e-6, vmin=1e-16, aspect='auto', origin="lower", interpolation='nearest', cmap="RdYlGn_r", extent=extent)
    
    im2 = axarr[1,i].imshow(np.sign(energyerror_sign[i]), vmax=1, vmin=-1, aspect='auto', origin="lower", interpolation='nearest', cmap="bwr", extent=extent)

    im3 = axarr[2,i].imshow(timing[i], vmin=0., vmax=3.*np.median(timing), aspect='auto', origin="lower", interpolation="nearest", cmap="RdYlGn_r", extent=extent)
        

try:  # Trips travis
    cax1,kw1 = matplotlib.colorbar.make_axes([ax for ax in axarr[0,:]])
    cax2,kw2 = matplotlib.colorbar.make_axes([ax for ax in axarr[1,:]])
    cax3,kw3 = matplotlib.colorbar.make_axes([ax for ax in axarr[2,:]])

    cb1 = plt.colorbar(im1, cax=cax1, **kw1)
    cb1.solids.set_rasterized(True)
    cb1.set_label("relative energy error")

    cb2 = plt.colorbar(im2, cax=cax2, **kw2)
    cb2.solids.set_rasterized(True)
    t = ticker.MaxNLocator(nbins=3)
    cb2.locator = t
    cb2.update_ticks()
    cb2.set_label("sign of energy error")

    cb3 = plt.colorbar(im3, cax=cax3, **kw3)
    cb3.solids.set_rasterized(True)
    cb3.set_label("runtime per timestep [$\mu$s]")
    cb3.update_ticks()
except:
    pass


plt.show()

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')
plt.savefig("2body.pdf",prop=fontP, bbox_inches='tight')
print "Average speedup (WH/Mikkola): %.4f" %(np.mean(timing[0])/np.mean(timing[1]))
from sys import platform as _platform
if _platform == "darwin":
    import os
    os.system("open 2body.pdf")
