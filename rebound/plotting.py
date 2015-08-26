# -*- coding: utf-8 -*-
import math
try:
    import matplotlib.pyplot as plt
    import numpy as np
    packages_missing = 0
except:
    packages_missing = 1
    pass

from .particle import Particle

def OrbitPlot(sim, figsize=(5,5), lim=None, Narc=100, unitlabel=None, color=False, showPeriastron=False, trails=False, lw=1.):
        if packages_missing == 1:
            printf("Matplotlib and/or numpy not found. Plotting functions not available\n")
            return
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        orbits = sim.calculate_orbits()
        particles = sim.particles
        x = [p.x for p in particles]
        y = [p.y for p in particles]
        if lim is None:
            lim = 0.
            for o in orbits:
                if o.a>0.:
                    r = (1.+o.e)*o.a
                else:
                    r = o.r
                if r>lim:
                    lim = r
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])

        if unitlabel is not None:
            unitlabel = " " + unitlabel
        else:
            unitlabel = ""
        ax.set_xlabel("x"+unitlabel)
        ax.set_ylabel("y"+unitlabel)

        if color:
            cm = plt.get_cmap("jet")
        else:
            cmf = plt.get_cmap("Greys")
            cm = lambda x: cmf(x/2.+0.5)


        ax.scatter(particles[0].x,particles[0].x, marker="*", s=35*lw, facecolor="black", edgecolor=None, zorder=3)
        for i, o in enumerate(orbits):
            primary = sim.calculate_com(i+1)
            colori = cm(float(i+1)/float(sim.N-1))
            pp = Particle(a=o.a, f=o.f, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
            ax.scatter(pp.x, pp.y, s=25*lw, facecolor="black", edgecolor=None, zorder=3)
            if o.a>0.: # bound orbit
                phase = np.linspace(0,2.*np.pi,Narc)
                for ph in phase:
                    newp = Particle(a=o.a, f=o.f+ph, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                    if trails:
                        alpha = ph/(2.*np.pi)
                        color = (colori[0], colori[1], colori[2], alpha)
                    else:
                        color = colori
                    ax.plot([pp.x, newp.x], [pp.y, newp.y], color=color, zorder=2, lw=lw)
                    pp = newp
            else:     # unbound orbit
                pp = None
                lim_phase = math.pi-math.acos(1./o.e)
                phase = np.linspace(-lim_phase*0.999,lim_phase*0.999,Narc)
                # add particle phase manually
                p_i = np.searchsorted(phase, o.f)
                phase = np.insert(phase, p_i, o.f)
                for ph in phase:
                    newp = Particle(a=o.a, f=ph, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                    if trails:
                        if ph<=o.f:
                            alpha = 1.-(-ph+o.f)/lim_phase
                        else:
                            alpha = 0.2
                        color = (colori[0], colori[1], colori[2], alpha)
                    else:
                        color = colori
                    if pp is not None:
                        ax.plot([pp.x, newp.x], [pp.y, newp.y], color=color, zorder=2, lw=lw)
                    pp = newp
            
            if showPeriastron:
                newp = Particle(a=o.a, f=0., inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                ax.plot([primary.x, newp.x], [primary.y, newp.y], linestyle="dotted", c=colori, zorder=1, lw=lw)
                ax.scatter([newp.x],[newp.y], marker="o", s=5.*lw, facecolor="none", edgecolor=colori, zorder=1)





        return fig


