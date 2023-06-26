from ctypes import c_uint32, c_uint, c_ulong, c_char_p, c_double, byref, CFUNCTYPE, POINTER
from . import clibrebound
import sys
import functools
import rebound

def hash(key):
    hash_types = c_uint32, c_uint, c_ulong
    PY3 = sys.version_info[0] == 3
    if PY3:
        string_types = str,
        int_types = int,
    else:
        string_types = basestring,
        int_types = int, long,
   
    if isinstance(key, int_types):
        return c_uint32(key)
    elif isinstance(key, hash_types):
        return key
    elif isinstance(key, string_types):
        clibrebound.reb_hash.restype = c_uint32
        return c_uint32(clibrebound.reb_hash(c_char_p(key.encode('ascii'))))
    else:
        raise AttributeError("Need to pash hash an integer or string.")

def mod2pi(x):
    try:
        x = float(x)
    except:
        ValueError("Argument of mod2pi needs to be a float.")
    clibrebound.reb_tools_mod2pi.restype = c_double
    return clibrebound.reb_tools_mod2pi(c_double(x))

def M_to_f(e, M):
    try:
        e = float(e)
        M = float(M)
    except:
        ValueError("Arguments of M_to_f need to be floats.")
    clibrebound.reb_tools_M_to_f.restype = c_double
    return clibrebound.reb_tools_M_to_f(c_double(e), c_double(M))

def E_to_f(e, E):
    try:
        e = float(e)
        E = float(E)
    except:
        ValueError("Arguments of E_to_f need to be floats.")
    clibrebound.reb_tools_E_to_f.restype = c_double
    return clibrebound.reb_tools_E_to_f(c_double(e), c_double(E))

def M_to_E(e, M):
    try:
        e = float(e)
        M = float(M)
    except:
        ValueError("Arguments of M_to_E need to be floats.")
    clibrebound.reb_tools_M_to_E.restype = c_double
    return clibrebound.reb_tools_M_to_E(c_double(e), c_double(M))

def spherical_to_xyz(magnitude=1., theta=0., phi=0.):
    """Initialize Cartesian vector from its magnitude and two spherical angles theta (polar angle measured from z) and phi (azimuthal angle measured from x)

    Arguments
    ---------
    magnitude: float
        Magnitude of the vector
    theta: float
        Polar angle of the vector, i.e. measured from the z axis
    phi: float
        Azimuthal angle of the vector, i.e. measured counterclockwise from the x axis

    Returns
    -------
    List of [x,y,z] components
    """ 
    clibrebound.reb_tools_spherical_to_xyz.restype = rebound._Vec3d
    xyz = clibrebound.reb_tools_spherical_to_xyz(c_double(magnitude), c_double(theta), c_double(phi))
    return rebound.Vec3d(xyz)

def xyz_to_spherical(vector):
    """Return spherical angle theta measured from z axis
    
    Arguments
    ---------
    vector: List-like (e.g., list, numpy array)
        3D Cartesian vector to convert to spherical coordinates

    Returns
    -------
    List [magnitude, theta, phi] of vector magnitude, polar angle theta (from z axis) and azimuthal angle phi (from x axis)
    """
    magnitude = c_double()
    theta = c_double()
    phi = c_double()
    clibrebound.reb_tools_xyz_to_spherical(rebound.Vec3d(vector)._vec3d, byref(magnitude), byref(theta), byref(phi))
    return magnitude.value, theta.value, phi.value

def deref_arg(contents_type, return_type=None):
    """
    Turns a Python function that accepts a value (e.g. Simulation) into a C
    function that accepts a pointer to a value.

    Used as a decorator factory.
    """
    _cfunctype = CFUNCTYPE(return_type, POINTER(contents_type))
    def _decorator(func):
        @functools.wraps(func)
        def _wrapper(pointer, *args, **kwargs):
            func(pointer.contents, *args, **kwargs)
        return _cfunctype(_wrapper)
    return _decorator
