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

def get_color(color):
    """
    Takes a string for a color name defined in matplotlib and returns of a 3-tuple of RGB values.
    Will simply return passed value if it's a tuple of length three.

    Parameters
    ----------
    color   : str
        Name of matplotlib color to calculate RGB values for.
    """

    if isinstance(color, tuple) and len(color) == 3: # already a tuple of RGB values
        return color

    try:
        import matplotlib.colors as mplcolors
    except:
        raise ImportError("Error importing matplotlib. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")
   
    try:
        hexcolor = mplcolors.cnames[color]
    except KeyError:
        raise AttributeError("Color not recognized in matplotlib.")

    hexcolor = hexcolor.lstrip('#')
    lv = len(hexcolor)
    return tuple(int(hexcolor[i:i + lv // 3], 16)/255. for i in range(0, lv, lv // 3)) # tuple of rgb values

def fading_line(x, y, color='black', alpha_initial=1., alpha_final=0., **kwargs):
    """
    Returns a matplotlib LineCollection connecting the points in the x and y lists, with a single color and alpha varying from alpha_initial to alpha_final along the line.
    Can pass any kwargs you can pass to LineCollection, like linewidgth.

    Parameters
    ----------
    x       : list or array of floats for the positions on the (plot's) x axis
    y       : list or array of floats for the positions on the (plot's) y axis
    color   : matplotlib color for the line. Can also pass a 3-tuple of RGB values (default: 'black')
    alpha_initial:  Limiting value of alpha to use at the beginning of the arrays.
    alpha_final:    Limiting value of alpha to use at the end of the arrays.
    """
    try:
        from matplotlib.collections import LineCollection
        from matplotlib.colors import LinearSegmentedColormap
        import numpy as np
    except:
        raise ImportError("Error importing matplotlib and/or numpy. Plotting functions not available. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")

    color = get_color(color)
    cdict = {'red': ((0.,color[0],color[0]),(1.,color[0],color[0])),
             'green': ((0.,color[1],color[1]),(1.,color[1],color[1])),
             'blue': ((0.,color[2],color[2]),(1.,color[2],color[2])),
             'alpha': ((0.,alpha_initial, alpha_initial), (1., alpha_final, alpha_final))}
    
    Npts = len(x)
    if len(y) != Npts:
        raise AttributeError("x and y must have same dimension.")
   
    segments = np.zeros((Npts-1,2,2))
    segments[0][0] = [x[0], y[0]]
    for i in range(1,Npts-1):
        pt = [x[i], y[i]]
        segments[i-1][1] = pt
        segments[i][0] = pt 
    segments[-1][1] = [x[-1], y[-1]]

    individual_cm = LinearSegmentedColormap('indv1', cdict)
    lc = LineCollection(segments, cmap=individual_cm, **kwargs)
    lc.set_array(np.linspace(0.,1.,len(segments)))
    return lc

def OrbitPlotOneSlice(sim, ax, lim=None, limz=None, Narc=100, color=False, periastron=False, trails=False, show_orbit=True, lw=1., axes="xy", plotparticles=[]):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap
    import numpy as np
    if not plotparticles: # default. No plotparticles list passed, so plot all obits
        ps = sim.particles[1:]
    else: 
        try:
            ps = [sim.particles[i] for i in plotparticles]
        except ParticleNotFound:
            raise AttributeError("Particle in plotparticles list passed to OrbitPlot not found in the simulation.")

    if lim is None:
        lim = 0.
        for p in ps:
            if p.a>0.:
                r = (1.+p.e)*p.a
            else:
                r = p.d
            if r>lim:
                lim = r
        lim *= 1.15
    if limz is None:
        z = [p.z for p in ps]
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
            colors = [get_color(color)]
        if isinstance(color, list):
            colors = []
            for c in color:
                colors.append(get_color(c))
    else:
        colors = ['black']
    coloriterator = cycle(colors)

    coords = {'x':0, 'y':1, 'z':2}
    axis0 = coords[axes[0]]
    axis1 = coords[axes[1]]
    
    ax.scatter(getattr(sim.particles[0],axes[0]),getattr(sim.particles[0],axes[1]), marker="*", s=35*lw, facecolor="black", edgecolor=None, zorder=3)
    
    proj = {}
    for p in ps:
        ax.scatter(getattr(p,axes[0]), getattr(p,axes[1]), s=25*lw, facecolor="black", edgecolor=None, zorder=3)
        colori = next(coloriterator)
       
        if show_orbit is True:
            alpha_final = 0. if trails is True else 1. # fade to 0 with trails

            hyperbolic = p.a < 0. # Boolean for whether orbit is hyperbolic
            if hyperbolic is False:
                o = np.array(p.sample_orbit())
                proj['x'],proj['y'],proj['z'] = [o[:,i] for i in range(3)]
                lc = fading_line(proj[axes[0]], proj[axes[1]], colori, alpha_final=alpha_final, lw=lw)
                ax.add_collection(lc)
            else:
                o = np.array(p.sample_orbit(useTrueAnomaly=False))
                # true anomaly stays close to limiting value and switches quickly at pericenter for hyperbolic orbit, so use mean anomaly
                proj['x'],proj['y'],proj['z'] = [o[:,i] for i in range(3)]
                lc = fading_line(proj[axes[0]], proj[axes[1]], colori, alpha_final=alpha_final, lw=lw)
                ax.add_collection(lc)
          
                alpha = 0.2 if trails is True else 1.
                o = np.array(p.sample_orbit(trailing=False, useTrueAnomaly=False))
                proj['x'],proj['y'],proj['z'] = [o[:,i] for i in range(3)]
                lc = fading_line(proj[axes[0]], proj[axes[1]], colori, alpha_initial=alpha, alpha_final=alpha, lw=lw)
                ax.add_collection(lc)

        if periastron:
            primary = p.jacobi_com
            newp = Particle(a=p.a, f=0., inc=p.inc, omega=p.omega, Omega=p.Omega, e=p.e, m=p.m, primary=primary, simulation=sim)
            ax.plot([getattr(primary,axes[0]), getattr(newp,axes[0])], [getattr(primary,axes[1]), getattr(newp,axes[1])], linestyle="dotted", c=colori, zorder=1, lw=lw)
            ax.scatter([getattr(newp,axes[0])],[getattr(newp,axes[1])], marker="o", s=5.*lw, facecolor="none", edgecolor=colori, zorder=1)

