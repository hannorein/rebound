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

def OrbitPlot(sim, figsize=(5,5), lim=None, Narc=100, unitlabel=None):
        if packages_missing == 1:
            printf("Matplotlib and/or numpy not found. Plotting functions not available\n")
            return
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        orbits = sim.calculate_orbits()
        particles = sim.particles
        x = [p.x for p in particles]
        y = [p.y for p in particles]
        r = [math.sqrt(p.x*p.x + p.y*p.y) for p in particles]
        if lim is None:
            lim = 1.4*max(r)
        ax.set_xlim([-lim,lim])
        ax.set_ylim([-lim,lim])

        if unitlabel is not None:
            unitlabel = " " + unitlabel
        else:
            unitlabel = ""
        ax.set_xlabel("x"+unitlabel)
        ax.set_ylabel("y"+unitlabel)

        xyo = np.zeros((2,Narc))
        phase = np.linspace(0,2.*np.pi,Narc)
        for i, o in enumerate(orbits):
            primary = sim.calculate_com(i+1)
            for j, ph in enumerate(phase):
                newp = Particle(a=o.a, anom=o.f+ph, MEAN=True, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                xyo[0,j] = newp.x
                xyo[1,j] = newp.y
            ax.plot(xyo[0], xyo[1])


        ax.scatter(x, y, s=25, c="lightgray")

        return fig


