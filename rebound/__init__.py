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

import sys
module = ReboundModule(__name__)
sys.modules[__name__] = module
