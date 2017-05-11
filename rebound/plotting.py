# -*- coding: utf-8 -*-
import math
from .particle import Particle
from itertools import cycle

def OrbitPlot(sim, figsize=None, lim=None, limz=None, Narc=100, unitlabel=None, color=False, periastron=False, trails=True, show_orbit=True, lw=1., slices=False, plotparticles=[]):
    """
    Convenience function for plotting instantaneous orbits.

    Parameters
    ----------
    slices          : bool, optional
        Plot all three slices if set to True. Default is False and plots orbits only in the xy plane.
    figsize         : tuple of float, optional
        Tuple defining the figure size (default: (5,5))
    lim             : float, optional           
        Limit for axes (default: None = automatically determined)
    limz            : float, optional           
        Limit for z axis, only used if slices=True (default: None = automatically determined)
    unitlabel       : str, optional          
        String describing the units, shown on axis labels (default: None)
    color           : bool, str or list, optional            
        By default plots in black. If set to True, plots using REBOUND color cycle. If a string or list of strings, e.g. ['red', 'cyan'], will cycle between passed colors.
    periastron  : bool, optional            
        Draw a marker at periastron (default: False)
    trails          : bool, optional            
        Draw trails instead of solid lines (default: False)
    show_orbit      : bool, optional
        Draw orbit trails/lines (default: True)
    lw              : float, optional           
        Linewidth (default: 1.)
    plotparticles   : list, optional
        List of particles to plot. Can be a list of any valid keys for accessing sim.particles, i.e., integer indices or hashes (default: plot all particles)

    Returns
    -------
    fig
        A matplotlib figure

    Examples
    --------
    The following example illustrates a typical use case.

    >>> sim = rebound.Simulation()
    >>> sim.add(m=1)
    >>> sim.add(a=1)
    >>> fig = rebound.OrbitPlot(sim)
    >>> fig.savefig("image.png") # save figure to file
    >>> fig.show() # show figure on screen

    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib import gridspec
        import numpy as np
    except:
        raise ImportError("Error importing matplotlib and/or numpy. Plotting functions not available. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")
    if unitlabel is not None:
        unitlabel = " " + unitlabel
    else:
        unitlabel = ""
    if slices:
        if figsize is None:
            figsize = (8,8)
        fig, ax = plt.subplots(2, 2, figsize=figsize)
        gs = gridspec.GridSpec(2, 2, width_ratios=[3., 2.], height_ratios=[2.,3.],wspace=0., hspace=0.) 
        OrbitPlotOneSlice(sim, plt.subplot(gs[2]), lim=lim, Narc=Narc, color=color, periastron=periastron, trails=trails, show_orbit=show_orbit, lw=lw, axes="xy", plotparticles=plotparticles)
        OrbitPlotOneSlice(sim, plt.subplot(gs[3]), lim=lim, limz=limz, Narc=Narc, color=color, periastron=periastron, trails=trails, show_orbit=show_orbit, lw=lw, axes="zy", plotparticles=plotparticles)
        OrbitPlotOneSlice(sim, plt.subplot(gs[0]), lim=lim, limz=limz, Narc=Narc, color=color, periastron=periastron, trails=trails, show_orbit=show_orbit, lw=lw, axes="xz", plotparticles=plotparticles)
        plt.subplot(gs[2]).set_xlabel("x"+unitlabel)
        plt.subplot(gs[2]).set_ylabel("y"+unitlabel)
      
        plt.setp(plt.subplot(gs[0]).get_xticklabels(), visible=False)
        plt.subplot(gs[0]).set_ylabel("z"+unitlabel)
        
        plt.subplot(gs[3]).set_xlabel("z"+unitlabel)
        plt.setp(plt.subplot(gs[3]).get_yticklabels(), visible=False)
    else:
        if figsize is None:
            figsize = (5,5)
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_xlabel("x"+unitlabel)
        ax.set_ylabel("y"+unitlabel)
        OrbitPlotOneSlice(sim, ax, lim=lim, Narc=Narc, color=color, periastron=periastron, trails=trails, show_orbit=show_orbit, lw=lw, plotparticles=plotparticles)
    return fig

def getcolor(color):
    import matplotlib.colors as mplcolors
    try:
        hexcolor = mplcolors.cnames[color]
    except KeyError:
        raise AttributeError("Color passed to OrbitPlot not recognized in matplotlib.")

    hexcolor = hexcolor.lstrip('#')
    lv = len(hexcolor)
    return tuple(int(hexcolor[i:i + lv // 3], 16)/255. for i in range(0, lv, lv // 3)) # tuple of rgb values

def OrbitPlotOneSlice(sim, ax, lim=None, limz=None, Narc=100, color=False, periastron=False, trails=False, show_orbit=True, lw=1., axes="xy", plotparticles=[]):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    if not plotparticles: # default. No plotparticles list passed, so plot all obits
        orbits = sim.calculate_orbits()
    else: 
        orbits = []
        for p in plotparticles:
            try:
                orbits.append(sim.particles[p].orbit)
            except ParticleNotFound:
                raise AttributeError("Particle in plotparticles list passed to OrbitPlot not found in the simulation.")

    particles = sim.particles
    if lim is None:
        lim = 0.
        for o in orbits:
            if o.a>0.:
                r = (1.+o.e)*o.a
            else:
                r = o.d
            if r>lim:
                lim = r
        lim *= 1.15
    if limz is None:
        z = [p.z for p in particles]
        limz = 2.0*max(z)
        if limz > lim:
            limz = lim
        if limz <= 0.:
            limz = lim


    if axes[0]=="z":
        ax.set_xlim([-limz,limz])
    else:
        ax.set_xlim([-lim,lim])
    if axes[1]=="z":
        ax.set_ylim([-limz,limz])
    else:
        ax.set_ylim([-lim,lim])

    if color:
        if color == True:
            colors = [(1.,0.,0.),(0.,0.75,0.75),(0.75,0.,0.75),(0.75, 0.75, 0,),(0., 0., 0.),(0., 0., 1.),(0., 0.5, 0.)]
        if isinstance(color, str):
            colors = [getcolor(color)]
        if isinstance(color, list):
            colors = []
            for c in color:
                colors.append(getcolor(c))
    else:
        colors = [(0.,0.,0.)]
    
    coloriterator = cycle(colors)

    ax.scatter(getattr(particles[0],axes[0]),getattr(particles[0],axes[1]), marker="*", s=35*lw, facecolor="black", edgecolor=None, zorder=3)
    for i, o in enumerate(orbits):
        colori = next(coloriterator)
        primary = sim.calculate_com(last=i)
        pp = Particle(a=o.a, f=o.f, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
        ax.scatter(getattr(pp,axes[0]), getattr(pp,axes[1]), s=25*lw, facecolor="black", edgecolor=None, zorder=3)
        if o.a>0.: # bound orbit
            segments = np.zeros((Narc,2,2))
            if show_orbit:
                phase = np.linspace(0,2.*np.pi,Narc)
                for j,ph in enumerate(phase):
                    newp = Particle(a=o.a, f=o.f+ph, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                    xy = [getattr(newp,axes[0]), getattr(newp,axes[1])]
                    segments[j][1]=xy
                    segments[(j+1)%Narc][0]=xy
            cdict = {'red': ((0.,colori[0],colori[0]),(1.,colori[0],colori[0])),
                     'green': ((0.,colori[1],colori[1]),(1.,colori[1],colori[1])),
                     'blue': ((0.,colori[2],colori[2]),(1.,colori[2],colori[2]))}
            if trails:
                cdict['alpha'] = ((0.,0.,0.),(1.,1.,1.))
            individual_cm = LinearSegmentedColormap('indv1', cdict)
            lc = LineCollection(segments, cmap=individual_cm, linewidth=lw)
            if show_orbit:
                lc.set_array(phase)
            ax.add_collection(lc)
        else:     # unbound orbit.  Step in M rather than f, since for hyperbolic orbits f stays near lim, and jumps to -f at peri
            t_cross = 4.*o.d/o.v # rough time to cross display box
            
            lim_phase = abs(o.n)*t_cross # M = nt, n is negative
            segments = np.zeros((Narc,2,2))
            if show_orbit:
                phase = np.linspace(-lim_phase+o.M,o.M,Narc)
                for j,ph in enumerate(phase):
                    newp = Particle(a=o.a, M=ph, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                    xy = [getattr(newp,axes[0]), getattr(newp,axes[1])]
                    segments[j][1]=xy
                    if j+1<Narc:
                       segments[j+1][0]=xy
                    if j==0:
                       segments[j][0]=xy
            cdict = {'red': ((0.,colori[0],colori[0]),(1.,colori[0],colori[0])),
                     'green': ((0.,colori[1],colori[1]),(1.,colori[1],colori[1])),
                     'blue': ((0.,colori[2],colori[2]),(1.,colori[2],colori[2]))}
            if trails:
                cdict['alpha'] = ((0.,0.,0.),(1.,1.,1.))
            individual_cm = LinearSegmentedColormap('indv1', cdict)
            lc = LineCollection(segments, cmap=individual_cm, linewidth=lw)
            if show_orbit:
                lc.set_array(phase)
            ax.add_collection(lc)
            
            segments = np.zeros((Narc,2,2))
            if show_orbit:
                phase = np.linspace(o.M,o.M+lim_phase,Narc)
                for j,ph in enumerate(phase):
                    newp = Particle(a=o.a, M=ph, inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
                    xy = [getattr(newp,axes[0]), getattr(newp,axes[1])]
                    segments[j][1]=xy
                    if j+1<Narc:
                       segments[j+1][0]=xy
                    if j==0:
                       segments[j][0]=xy
            cdict = {'red': ((0.,colori[0],colori[0]),(1.,colori[0],colori[0])),
                     'green': ((0.,colori[1],colori[1]),(1.,colori[1],colori[1])),
                     'blue': ((0.,colori[2],colori[2]),(1.,colori[2],colori[2]))}
            if trails:
                cdict['alpha'] = ((0.,0.2,0.2),(1.,0.2,0.2))
            individual_cm = LinearSegmentedColormap('indv1', cdict)
            lc = LineCollection(segments, cmap=individual_cm, linewidth=lw)
            if show_orbit:
                lc.set_array(phase)
            ax.add_collection(lc)
        
        
        if periastron:
            newp = Particle(a=o.a, f=0., inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=particles[i+1].m, primary=primary, simulation=sim)
            ax.plot([getattr(primary,axes[0]), getattr(newp,axes[0])], [getattr(primary,axes[1]), getattr(newp,axes[1])], linestyle="dotted", c=colori, zorder=1, lw=lw)
            ax.scatter([getattr(newp,axes[0])],[getattr(newp,axes[1])], marker="o", s=5.*lw, facecolor="none", edgecolor=colori, zorder=1)

