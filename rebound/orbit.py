import ctypes
from .tools import M_to_E
from .vectors import Vec3dBasic

class Orbit(ctypes.Structure):
    """
    A class containing orbital parameters for a particle.
    This is an abstraction of the reb_orbit data structure in C.

    When using the various REBOUND functions using Orbits, all angles are in radians. 
    The following image illustrated the most important angles used.
    In REBOUND the reference direction is the positive x direction, the reference plane
    is the xy plane.
    
    .. image:: images/orbit.png
       :width: 500px
       :height: 450px

    Image from wikipedia. CC-BY-SA-3.

    Attributes
    ----------
    d       : float           
        radial distance from reference 
    v       : float         
        velocity relative to central object's velocity
    h       : float           
        specific angular momentum
    P       : float           
        orbital period (negative if hyperbolic)
    n       : float           
        mean motion    (negative if hyperbolic)
    a       : float           
        semimajor axis
    e       : float           
        eccentricity
    inc     : float           
        inclination
    Omega   : float           
        longitude of ascending node
    omega   : float           
        argument of pericenter
    pomega  : float           
        longitude of pericenter
    f       : float           
        true anomaly
    M       : float           
        mean anomaly
    E       : float           
        eccentric anomaly (requires solving Kepler's equation - only calculated when needed)
    l       : float           
        mean longitude = Omega + omega + M
    theta   : float           
        true longitude = Omega + omega + f
    T       : float
        time of pericenter passage
    rhill   : float
        Hill radius ( =a*pow(m/(3M),1./3.) )
    pal_h   : float
        Cartesian component of the eccentricity ( = e*sin(pomega) )
    pal_k   : float
        Cartesian component of the eccentricity ( = e*cos(pomega) )
    pal_ix   : float
        Cartesian component of the inclination ( = 2*sin(i/2)*cos(Omega) )
    pal_iy   : float
        Cartesian component of the inclination ( = 2*sin(i/2)*sin(Omega) )
    hvec     : Vec3d
        Specific angular momentum vector
    evec     : Vec3d
        Eccentricity vector (mag = eccentricity, points toward pericenter)
    """
    _fields_ = [("d", ctypes.c_double),
                ("v", ctypes.c_double),
                ("h", ctypes.c_double),
                ("P", ctypes.c_double),
                ("n", ctypes.c_double),
                ("a", ctypes.c_double),
                ("e", ctypes.c_double),
                ("inc", ctypes.c_double),
                ("Omega", ctypes.c_double),
                ("omega", ctypes.c_double),
                ("pomega", ctypes.c_double),
                ("f", ctypes.c_double),
                ("M", ctypes.c_double),
                ("l", ctypes.c_double),
                ("theta", ctypes.c_double),
                ("T", ctypes.c_double),
                ("rhill", ctypes.c_double),
                ("pal_h", ctypes.c_double),
                ("pal_k", ctypes.c_double),
                ("pal_ix", ctypes.c_double),
                ("pal_iy", ctypes.c_double),
                ("hvec", Vec3dBasic),
                ("evec", Vec3dBasic)]

    def __repr__(self):
        return "<rebound.Orbit instance, a={0} e={1} inc={2} Omega={3} omega={4} f={5}>".format(str(self.a),str(self.e), str(self.inc), str(self.Omega), str(self.omega), str(self.f))

    @property 
    def E(self):
        return M_to_E(self.e, self.M)

