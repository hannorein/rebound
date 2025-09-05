import ctypes
import math
from .tools import M_to_E
from .vectors import Vec3dBasic
from .tools import E_to_f, M_to_f, mod2pi
from .particle import Particle

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
    
    def sample(self, primary, Npts=100, samplingAngle=None, duplicateEndpoint=None):
        """
        Returns a nested list of xyz positions along the osculating orbit of the particle. 
        If primary is not passed, returns xyz positions along the Jacobi osculating orbit
        (with mu = G*Minc, where Minc is the total mass from index 0 to the particle's index, inclusive). 

        Parameters
        ----------
        Npts    : int, optional  
            Number of points along the orbit to return  (default: 100)
        primary : rebound.Particle, optional
            Primary to use for the osculating orbit (default: Jacobi center of mass)
        samplingAngle: str, optional
            This determines which angle is sampled linearly. Can be:
              - "M" (mean anomaly)
              - "E" (eccentric anomaly) 
              - "f" (true anomaly)
        duplicateEndpoint: bool, optional
            If true, then the first and last point will be identical for closed orbits. This is useful for some plotting tools. The default is true for eccentric orbits. The argument has no effect for hyperbolic orbits (because the endpoints are not identical).
        """

        phases_f = []
        if samplingAngle is not None:
            if any(c not in "EMf" for c in samplingAngle):
                raise ValueError("Unknown character in samplingAngle.")
        
        if self.a < 0.: # hyperbolic orbit
            if samplingAngle is None:
                samplingAngle = "Mf"
            Nptsangle = {}
            for angle in samplingAngle[1:]:
                Nptsangle[angle] = (Npts-1)//len(samplingAngle) # one point is reserved for actual position
            Nptsangle[samplingAngle[0]] = Npts-1-sum(Nptsangle.values())
            if "M" in samplingAngle:
                phi = math.acos(-1./self.e)*0.999
                Npts = Nptsangle["M"]
                dphi = 2*phi/(Npts-1)
                for i in range(Npts):
                    f = M_to_f(self.e, phi)
                    phases_f.append(f)
                    phi -= dphi
            if "E" in samplingAngle:
                phi = math.acos(-1./self.e)*0.999
                Npts = Nptsangle["E"]
                dphi = 2*phi/(Npts-1)
                for i in range(Npts):
                    f = E_to_f(self.e, phi)
                    phases_f.append(f)
                    phi -= dphi
            if "f" in samplingAngle:
                phi = math.acos(-1./self.e)*0.999
                Npts = Nptsangle["f"]
                dphi = 2*phi/(Npts-1)
                for i in range(Npts):
                    f = mod2pi(phi)
                    phases_f.append(f)
                    phi -= dphi
        else:       # circular orbit
            if samplingAngle is None:
                samplingAngle = "Ef"
            if duplicateEndpoint is None:
                duplicateEndpoint = True
            Nptsangle = {}
            for angle in samplingAngle[1:]:
                Nptsangle[angle] = (Npts-1)//len(samplingAngle) # one point is reserved for actual position
            Nptsangle[samplingAngle[0]] = Npts-1-sum(Nptsangle.values())
            if "M" in samplingAngle:
                Npts = Nptsangle["M"]
                dphi = 2.*math.pi/(Npts-1 if duplicateEndpoint else Npts)  # one point is reserved for the end point
                for i in range(Npts):
                    f = M_to_f(self.e, i*dphi)
                    phases_f.append(f)
            if "E" in samplingAngle:
                Npts = Nptsangle["E"]
                dphi = 2.*math.pi/(Npts-1 if duplicateEndpoint else Npts)  # one point is reserved for the end point
                for i in range(Npts):
                    f = E_to_f(self.e, i*dphi)
                    phases_f.append(f)
            if "f" in samplingAngle:
                Npts = Nptsangle["f"]
                dphi = 2.*math.pi/(Npts-1 if duplicateEndpoint else Npts)  # one point is reserved for the end point
                for i in range(Npts):
                    f = i*dphi
                    f = mod2pi(f)
                    phases_f.append(f)

        # add actual position
        f = mod2pi(self.f)
        phases_f.append(f)
        phases_f.sort()
      
        pts_pre = []
        pts_post = []
        
        for f in phases_f:
            newp = Particle(a=self.a, f=f, inc=self.inc, omega=self.omega, Omega=self.Omega, e=self.e, primary=primary, simulation=primary._sim.contents)
            if f<=self.f:
                pts_pre.append(newp.xyz)
            else:
                pts_post.append(newp.xyz)
        
        return pts_post + pts_pre

