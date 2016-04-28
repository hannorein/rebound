# -*- coding: utf-8 -*-

"""
Initial conditions for standard tests

"""

import math

def add_outer_solar_system(sim):
    """
    Add the planet of the outer Solar System as a test problem.
    Data taken from NASA Horizons.
    """
    Gfac = 1./0.01720209895       # Gaussian constant 
    if sim.G is not None:
        Gfac *= math.sqrt(sim.G)

    sim.add( m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6*Gfac, vy=-6.33922479583593e-6*Gfac, vz=-3.13202145590767e-9*Gfac )  # Sun
    sim.add( m=1./1047.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3*Gfac, vy=+5.51815399480116e-3*Gfac, vz=-2.66711392865591e-6*Gfac )  # Jupiter
    sim.add( m=1./3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3*Gfac, vy=+3.99723751748116e-3*Gfac, vz=+1.67206320571441e-5*Gfac )  # Saturn
    sim.add( m=1./22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3*Gfac, vy=+2.06438412905916e-3*Gfac, vz=-2.17699042180559e-5*Gfac )  # Uranus
    sim.add( m=1./19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4*Gfac, vy=-3.11361111025884e-3*Gfac, vz=+3.58344705491441e-5*Gfac )  # Neptune
    sim.add( m=7.4074074e-09, x=-2.13858977531573e+1, y=+3.20719104739886e+1, z=+2.49245689556096e+0, vx=-1.76936577252484e-3*Gfac, vy=-2.06720938381724e-3*Gfac, vz=+6.58091931493844e-4*Gfac )  # Pluto


