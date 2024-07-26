# -*- coding: utf-8 -*-
from .particle import Particle
from itertools import cycle


class OrbitPlot:
    """
    Class for visualizing simulations using instantaneous orbits.
    """

    def __init__(self, sim, fig=None, ax=None, figsize=(5,5), projection="xy", xlim=None, ylim=None, unitlabel=None, color=False, periastron=False, orbit_style="trail", lw=1., particles=None, primary=None, show_primary=True, origin=None, Narc=128):
        """
        Initializer for OrbitPlot class

        By default each instance of OrbitPlot creates its own matplotlib figure and axes. However, you can also reuse an existing figure and axes by passing them as arguments. 

        Parameters
        ----------
        sim             : Simulation (required)
        fig             : matplotlib.figure.Figure, optional
            If a figure instances is passed as an argument, we will use it instead of creating a new one.
        ax              : matplotlib.axes._subplots.AxesSubplot, optional
            If a Axes instances is passed as an argument, we will use it instead of creating a new one.
        figsize         : tuple of float, optional
            Tuple defining the figure size (default: (5,5))

        projection      : string, optional
            By default the orbit is shown as a projection into the xy plane. To show the orbit projected into the xz plane, set this string to "xy".
        xlim            : tuple of float, optional           
            Limits for x axes (default: None = automatically determined)
        ylim            : tuple of float, optional           
            Limits for y axes (default: None = automatically determined)
        unitlabel       : str, optional          
            String describing the units, shown on axis labels (default: None)

        color           : bool, str or list, optional            
            By default plots are black and white. If set to True, plots use a color cycle. If set to a string or list of strings, e.g. ['red', 'cyan'], OrbitPlot will cycle between the colors.
        periastron  : bool, optional            
            Draw a marker at periastron (default: False)
        orbit_style     : str, optional
            This argument determines the type of orbit show. By default, it shows the orbit as a trailing and fading line ("trail"). Other object are: "solid", None.
        lw              : float, optional           
            Linewidth used in plots (default: 1.)

        particles       : list of (int, or str), optional
            List of particles to plot. Can be a list of any valid keys for accessing sim.particles, i.e., integer indices or hashes (default: plot all particles). List should not include primary.
        primary         : rebound.Particle, optional
            Primary to use for the osculating orbit (default: Jacobi center of mass)
        show_primary    : bool
            Set to False to hide primary (default: True)
        origin          : tuple of two floats, optional
            By default the origin to [0, 0]. Use this argument if you wish to move the entire OrbitPlot.
        Narc            : int, optional
            Number of points used in an orbit. Increase this number for highly eccentric orbits. (default: 128)

        Returns
        -------
        The function returns a new instance of the OrbitPlot class. 

        Examples
        --------
        The following example illustrates a typical use case.

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1)
        >>> sim.add(a=1)
        >>> op = rebound.OrbitPlot(sim)
        >>> op.fig.savefig("image.png") # save figure to file

        """
        try:
            import matplotlib.pyplot as plt
        except:
            raise ImportError("Error importing matplotlib and/or numpy. Plotting functions not available. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")

        self.sim = sim 
        
        self._xlim = xlim
        self._ylim = ylim
        self.color = color
        self.projection = projection
        self._periastron = periastron # cannot be changed later
        self._orbit_style = orbit_style # cannot be changed later
        self._lw = lw # cannot be changed later
        if isinstance(primary,(str,int)):
            primary = sim.particles[primary]
        self._primary = primary # cannot be changed later
        self._show_primary = show_primary # cannot be changed later
        self._Narc = Narc # cannot be changed later
        self._particles = particles # cannot be changed later
        self.origin = origin
           
        updateLimits = True
        if fig is not None:
            updateLimits = False
            self.fig = fig
            self.ax = ax
        else:
            # Initialize figure (figure size can not be changed later)
            if unitlabel is not None:
                unitlabel = " " + unitlabel
            else:
                unitlabel = ""
            self.fig = plt.figure(figsize=figsize)
            self.ax = plt.subplot(111,aspect="equal") #check syntax
            self.ax.set_xlabel("x"+unitlabel)
            self.ax.set_ylabel("y"+unitlabel)

        self.orbits = None
        self.primary = None
        self.particles = None
        self.draw(updateLimits=updateLimits)

    @property
    def xlim(self):
        if self._xlim:
            return self._xlim
        else:
            return self.ax.get_xlim()
    
    @xlim.setter
    def xlim(self, value):
        self._xlim = value
    
    @property
    def ylim(self):
        if self._ylim:
            return self._ylim
        else:
            return self.ax.get_ylim()
    
    @ylim.setter
    def ylim(self, value):
        self._ylim = value


    def draw(self, update=False, updateLimits=True):
        if self.particles is None or update==False:
            update = True # First run needs update
            self.setup()
        if update:
            self.update(updateLimits=updateLimits)

    def offset(self):
        import numpy as np
        if self.origin is not None:
            if isinstance(self.origin, (list, np.ndarray)):
                return [-self.origin[0],-self.origin[1]]
            elif isinstance(self.origin, Particle):
                p = self.origin    
            else:
                p = self.sim.particles[self.origin]
            px, py = getattr(p,self.projection[0]), getattr(p,self.projection[1])
            return [-px, -py]
        else:
            return [0, 0]

    def setup(self):
        """
        Construct the collections for particles and orbits. 
        Does not populate the data arrays -- that is done in update().
        """
        from matplotlib.collections import LineCollection
        import numpy as np

        particles = self._particles
        if not particles:
            particles = range(1, self.sim.N_real)

        # Color stuff
        color = self.color
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
             colors = [(0.,0.,0.)]
        coloriterator = cycle(colors)
        color_list = []
        for i in range(len(particles)):
            colori = next(coloriterator)
            color_list.append(colori)

        if self._show_primary:
            pc = self.ax.scatter([],[], marker="*", s=35*self._lw, facecolor="black", edgecolor=None, zorder=3)
            self.ax.add_collection(pc)
            self.primary = pc
        
        pc = self.ax.scatter([], [], s=25*self._lw, facecolor="black", edgecolor=None, zorder=3)
        self.ax.add_collection(pc)
        self.particles = pc
        
        lcs = []

        for j in range(len(particles)):
            if self._orbit_style is not None:
                lc = LineCollection([], lw=self._lw)

                line_colors = np.zeros((self._Narc+1,4))
                line_colors[:,0:3] = color_list[j]
                alpha = 1.
                if self._orbit_style=="trail":
                    line_colors[:,3] = alpha*np.linspace(0,1,self._Narc+1)
                elif self._orbit_style=="solid":
                    line_colors[:,3] = alpha
                else:
                    raise ValueError("Unknown orbit_style.")
                lc.set_color(line_colors)

                lcs.append(lc)
                self.ax.add_collection(lc)

        self.orbits = lcs

        if self._periastron:
            lc = LineCollection([], lw=self._lw, zorder=1, linestyle="dotted")
            lc.set_color(color_list)
            self.ax.add_collection(lc)
            self.periastrons = lc


    def update(self, updateLimits=False):
        """
        This function sets the data arrays in the plot according to the simulation in self.sim.
        Needs to be called after setup() to actually show something on the figure.
        """

        if self._xlim and self._ylim:
            updateLimits = False # No need to calculate limits.

        offset_x, offset_y = self.offset()
        
        import numpy as np
        projection = self.projection
        particles = self._particles
        if not particles:
            particles = range(1, self.sim.N_real)
        
        prim = self.sim.particles[0] if self._primary is None else self._primary 
        px, py = getattr(prim,projection[0])+offset_x, getattr(prim,projection[1])+offset_y
        if self._show_primary:
            pc = self.primary
            pc.set_offsets([[px, py]])
        limits_x = [px, px]
        limits_y = [py, py]

        offsets = np.zeros((len(particles),2))
        periastrons = []
        proj = {"x": 0, "y": 1, "z": 2}
        for j in range(len(particles)):
            p = self.sim.particles[particles[j]]
            px, py = getattr(p,projection[0])+offset_x, getattr(p,projection[1])+offset_y
            offsets[j] = [px, py]
            
            if len(self.orbits): # no need to update if orbit_style = None
                # Need many line segments here so we can color it.
                pts = np.array(p.sample_orbit(Npts=self._Narc+1, primary=self._primary))
                segments = np.zeros((self._Narc,2,2))
                x, y = pts[:,proj[projection[0]]], pts[:,proj[projection[1]]]
                segments[:,0,0] = x[:-1] + offset_x
                segments[:,0,1] = y[:-1] + offset_y
                segments[:,1,0] = x[1:] + offset_x
                segments[:,1,1] = y[1:] + offset_y
                lc = self.orbits[j]
                lc.set_segments(segments)
           
                if updateLimits:
                    if p.a < 0.: # hyperbolic, avoid zooming out too much
                        limits_x = [min(limits_x[0], px), max(limits_x[1], px)] 
                        limits_y = [min(limits_y[0], py), max(limits_y[1], py)] 
                    else:
                        ma = np.max(segments, axis=(0,1))
                        mi = np.min(segments, axis=(0,1))
                        limits_x = [min(limits_x[0], px, mi[0]), max(limits_x[1], px, ma[0])] 
                        limits_y = [min(limits_y[0], py, mi[1]), max(limits_y[1], py, ma[1])] 
                    

                if self._periastron:
                    if self._primary is None:
                        pprimary = p.jacobi_com
                    else:
                        pprimary = self._primary
                    o = p.orbit(primary = pprimary)
                    newp = Particle(a=o.a, f=0., inc=o.inc, omega=o.omega, Omega=o.Omega, e=o.e, m=p.m, primary=pprimary, simulation=self.sim)
                    periastrons.append([[getattr(prim,projection[0])+offset_x, getattr(prim,projection[1])+offset_y], 
                                        [getattr(newp,projection[0])+offset_x, getattr(newp,projection[1])+offset_y]])
            else:
                if updateLimits:
                    limits_x = [min(limits_x[0], px), max(limits_x[1], px)] 
                    limits_y = [min(limits_y[0], py), max(limits_y[1], py)] 

        self.particles.set_offsets(offsets) # path collection offsets
        if self._periastron:
            self.periastrons.set_segments(periastrons)
            
        if updateLimits:
            width  = limits_x[1] - limits_x[0]
            height = limits_y[1] - limits_y[0]
        
            # prevent overly elongated plots
            limits_x[0] -= height/10.
            limits_x[1] += height/10.
            limits_y[0] -= width/10.
            limits_y[1] += width/10.

            if self._xlim is None:
                self.ax.set_xlim([limits_x[0]-0.05*width, limits_x[1]+0.05*width ])
            if self._ylim is None:
                self.ax.set_ylim([limits_y[0]-0.05*height,limits_y[1]+0.05*height])

        if self._xlim:
            self.ax.set_xlim(self._xlim)
        if self._ylim:
            self.ax.set_ylim(self._ylim)
        

class OrbitPlotSet:
    """
    Class for visualizing simulations using instantaneous orbits in 3D. Uses three rebound.OrbitPlot instances internally.
    """

    def __init__(self, sim, slices=0.5, fig=None, ax=None, figsize=(5,8), unitlabel=None, **kwargs):
        """
        Initializer for OrbitPlotSet class

        This function has the same arguments as OrbitPlot() with the addition of:

        Parameters
        ----------
        slices            : float, optional
            Changes the height and width of the top and right plot relative to the main plot. Default: 0.5.
        """

        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.axes_grid1 import make_axes_locatable
        except:
            raise ImportError("Error importing matplotlib and/or numpy. Plotting functions not available. If running from within a jupyter notebook, try calling '%matplotlib inline' beforehand.")
        updateLimits = True
        if fig is not None:
            updateLimits = False
            self.fig = fig
            self.ax_main, self.ax_top, self.ax_right = ax
        else:
            # Initialize figure (figure size can not be changed later)
            
            if unitlabel is not None:
                unitlabel = " " + unitlabel
            else:
                unitlabel = ""

            self.fig = plt.figure(figsize=figsize)
            self.ax_main = plt.subplot(111,aspect="equal") #check syntax
            self.ax_main.set_xlabel("x"+unitlabel)
            self.ax_main.set_ylabel("y"+unitlabel)

            divider = make_axes_locatable(self.ax_main)
            divider.set_aspect(True)
            self.ax_top   = divider.append_axes("top",  size="%.2f%%"%(100.*slices), sharex=self.ax_main, pad=0)
            self.ax_top.set_aspect('equal', adjustable='datalim')
            self.ax_right = divider.append_axes("right", size="%.2f%%"%(100.*slices), sharey=self.ax_main, pad=0)
            self.ax_right.set_aspect('equal', adjustable='datalim')
          
            plt.setp(self.ax_top.get_xticklabels(), visible=False)
            plt.setp(self.ax_top.get_xticklines(), visible=False)
            self.ax_top.set_ylabel("z"+unitlabel)
            
            plt.setp(self.ax_right.get_yticklabels(), visible=False)
            plt.setp(self.ax_right.get_yticklines(), visible=False)
            self.ax_right.set_xlabel("z"+unitlabel)

        self.sim = sim
        self.main  = OrbitPlot(sim, fig=self.fig, ax=self.ax_main, projection="xy", **kwargs)
        self.top   = OrbitPlot(sim, fig=self.fig, ax=self.ax_top, projection="xz", **kwargs)
        self.right = OrbitPlot(sim, fig=self.fig, ax=self.ax_right, projection="zy", **kwargs)
        
        self.draw(updateLimits=updateLimits, update=True)

    def draw(self, update=False, updateLimits=True):
        self.main.draw(update=update, updateLimits=updateLimits)
        self.top.draw(update=update, updateLimits=updateLimits)
        self.right.draw(update=update, updateLimits=updateLimits)
    def update(self, updateLimits=True):
        self.main.update(updateLimits=updateLimits)
        self.top.update(updateLimits=updateLimits)
        self.right.update(updateLimits=updateLimits)


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

    import matplotlib.colors as mplcolors
   
    try:
        hexcolor = mplcolors.cnames[color]
    except KeyError:
        raise AttributeError("Color not recognized in matplotlib.")

    hexcolor = hexcolor.lstrip('#')
    lv = len(hexcolor)
    return tuple(int(hexcolor[i:i + lv // 3], 16)/255. for i in range(0, lv, lv // 3)) # tuple of rgb values

#def OrbitPlotAddFancyStars(ax,lw,slices=1.):
#    import numpy as np
#    # safe the current random seed to restore later
#    os = np.random.get_state()
#    # always produce the same stars
#    np.random.seed(1) 
#
#    ax.set_facecolor((0.,0.,0.))
#    for pos in ['top', 'bottom', 'right', 'left']:
#        ax.spines[pos].set_edgecolor((0.3,0.3,0.3))
#    
#    starcolor = (1.,1.,1.)
#    starsurfacedensity = 0.8
#
#    area = np.sqrt(np.sum(np.square(ax.transAxes.transform([1.,1.]) - ax.transAxes.transform([0.,0.]))))*slices
#    nstars = int(starsurfacedensity*area)
#
#    #small stars
#    xy = np.random.uniform(size=(nstars,2))
#    ax.scatter(xy[:,0],xy[:,1], transform=ax.transAxes, alpha=0.05, s=8*lw, facecolor=starcolor, edgecolor=None, zorder=3)
#    ax.scatter(xy[:,0],xy[:,1], transform=ax.transAxes, alpha=0.1, s=4*lw, facecolor=starcolor, edgecolor=None, zorder=3)
#    ax.scatter(xy[:,0],xy[:,1], transform=ax.transAxes, alpha=0.2, s=0.5*lw, facecolor=starcolor, edgecolor=None, zorder=3)
#    
#    #large stars
#    xy = np.random.uniform(size=(nstars//4,2))
#    ax.scatter(xy[:,0],xy[:,1], transform=ax.transAxes, alpha=0.1, s=15*lw, facecolor=starcolor, edgecolor=None, zorder=3)
#    ax.scatter(xy[:,0],xy[:,1], transform=ax.transAxes, alpha=0.1, s=5*lw, facecolor=starcolor, edgecolor=None, zorder=3)
#    ax.scatter(xy[:,0],xy[:,1], transform=ax.transAxes, alpha=0.5, s=2*lw, facecolor=starcolor, edgecolor=None, zorder=3)
#
#    np.random.set_state(os)
#
#
#
#      
#        if self.fancy:
#            colors = [(181./206.,66./206.,191./206.)]
#        else:
#        if self.fancy:
#            #TODO Color=colori
#            pc = self.ax.scatter([],[], s=25*self._lw, edgecolor=None, zorder=3) # need to set color manually later in update
##        if self.fancy:
#            prim = self.sim.particles[0] if self._primary is None else self._primary 
#            #TODO:
#            sun = (256./256.,256./256.,190./256.)
#            opa = 0.035
#            size = 6000.
#            for i in range(100):
#                ax.scatter(getattr(prim,axes[0]),getattr(prim,axes[1]), alpha=opa, s=size*self._lw, facecolor=sun, edgecolor=None, zorder=3)
#                size *= 0.95
#            pc = self.ax.scatter(getattr(prim,axes[0]),getattr(prim,axes[1]), s=size*self._lw, facecolor=sun, edgecolor=None, zorder=3)
#        else:
