# -*- coding: utf-8 -*-

"""
REBOUND Module. 

"""
#Make changes for python 2 and 3 compatibility
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x


from .librebound import ReboundModule
from .particle import Particle
from .particle import Orbit
from .interruptible_pool import InterruptiblePool
import sys
module = ReboundModule(__name__)
module.Particle = Particle
module.Orbit = Orbit
module.InterruptiblePool = InterruptiblePool
sys.modules[__name__] = module
