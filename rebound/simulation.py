from ctypes import Structure, c_double, POINTER, c_uint32, c_float, c_int, c_uint, c_uint32, c_int64, c_long, c_ulong, c_ulonglong, c_void_p, c_char_p, CFUNCTYPE, byref, create_string_buffer, addressof, pointer, cast, c_char, c_size_t, string_at
from . import clibrebound, Escape, NoParticles, Encounter, Collision, SimulationError, ParticleNotFound, M_to_E
from .citations import cite
from .particle import Particle
from .units import units_convert_particle, check_units, convert_G, hash_to_unit
from .tools import hash as rebhash
import math
import os
import sys
import ctypes.util
import warnings
try:
    # Required for Python>=3.9
    from collections.abc import MutableMapping
except:
    from collections import MutableMapping


import types
      
### The following enum and class definitions need to
### consitent with those in rebound.h
        
INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "leapfrog": 4, "none": 7, "janus": 8, "mercurius": 9, "saba": 10, "eos": 11, "bs": 12, "whfast512":21}
BOUNDARIES = {"none": 0, "open": 1, "periodic": 2, "shear": 3}
GRAVITIES = {"none": 0, "basic": 1, "compensated": 2, "tree": 3, "mercurius": 4, "jacobi": 5}
COLLISIONS = {"none": 0, "direct": 1, "tree": 2, "mercurius": 3, "line": 4, "linetree": 5}
VISUALIZATIONS = {"none": 0, "opengl": 1, "webgl": 2}
WHFAST_KERNELS = {"default": 0, "modifiedkick": 1, "composition": 2, "lazy": 3}
WHFAST_COORDINATES = {"jacobi": 0, "democraticheliocentric": 1, "whds": 2}
SABA_TYPES = {
        "1": 0x0, "2": 0x1, "3": 0x2, "4": 0x3,
        "cm1": 0x100, "cm2": 0x101, "cm3": 0x102, "cm4": 0x103,
        "cl1": 0x200, "cl2": 0x201, "cl3": 0x202, "cl4": 0x203,
        "10,4": 0x4, "8,6,4": 0x5, "10,6,4": 0x6,
        "h8,4,4": 0x7, "h8,6,4": 0x8, "h10,6,4": 0x9,
        }
EOS_TYPES = {
        "lf": 0x00,
        "lf4": 0x01,
        "lf6": 0x02,
        "lf8": 0x03,
        "lf4_2": 0x04,
        "lf8_6_4": 0x05,
        "plf7_6_4": 0x06,
        "pmlf4": 0x07,
        "pmlf6": 0x08,
        }

# Format: Majorerror, id, message
BINARY_WARNINGS = [
    (True,  1, "Cannot read binary file. Check filename and file contents."),
    (False, 2, "Binary file was saved with a different version of REBOUND. Binary format might have changed."),
    (False, 4, "You have to reset function pointers after creating a reb_simulation struct with a binary file."),
    (False, 8, "Binary file might be corrupted. Number of particles found does not match particle number expected."),
    (True,  16, "Error while reading binary file (file was closed).",),
    (True,  32, "Index out of range.",),
    (True,  64, "Error while trying to seek file.",),
    (False, 128, "Encountered unkown field in file. File might have been saved with a different version of REBOUND."),
    (True,  256, "Integrator type is not supported by this simulation archive version."),
    (False,  512, "The binary file seems to be corrupted. An attempt has been made to read the uncorrupted parts of it."),
    (True, 1024, "Reading old SimulationArchives (version < 2) is no longer supported. If you need to read such an archive, use a REBOUND version <= 3.26.3"),
]

class reb_hash_pointer_pair(Structure):
    _fields_ = [("hash", c_uint32),
                ("index", c_int)]

class reb_binary_field_descriptor(Structure):
    def __repr__(self):
        return '<{0}.{1} object at {2}, type={3}, dtype={4}, name=\'{5}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.type, self.dtype, self.name.decode("ascii"))
    _fields_ = [("type", c_uint),
                ("dtype", c_int),
                ("name", c_char*1024),
                ("offset", c_size_t),
                ("offset_N", c_size_t),
                ("element_size", c_size_t),
                ]
def binary_field_descriptor_list():
    fd_pointer_t = POINTER(reb_binary_field_descriptor)
    fd_pointer = (reb_binary_field_descriptor*3).in_dll(clibrebound, "reb_binary_field_descriptor_list")
    fd_pointer = cast(fd_pointer, fd_pointer_t) # not sure why I have to do it this way
    l = []
    i=0
    while True:
        l.append(fd_pointer[i])
        if fd_pointer[i].name == b'end':
            break
        i += 1
    return l


class Rotation(Structure):
    """
    This class facilitates rotations of Vec3d objects, and provides various convenience functions
    for commonly used rotations in celestial mechanics.
    """
    def __init__(self, ix=None, iy=None, iz=None, r=None, angle=None, axis=None):
        """
        Rotations are implemented as quaternions r + (ix)i + (iy)j + (iz)k. To initialize one
        can directly pass a set of the real numbers (ix, iy, iz, r). Alternatively one can pass
        an axis vector for the rotation axis, and an angle of rotation (counter-clockwise around axis).
        Only one full set of (r, ix, iy, iz) OR (angle, axis) must be passed.

        Arguments
        ---------
        ix : float
            Coefficient of i in quaternion r + (ix)i + (iy)j + (iz)k
        iy : float
            Coefficient of j in quaternion r + (ix)i + (iy)j + (iz)k
        iz : float
            Coefficient of k in quaternion r + (ix)i + (iy)j + (iz)k
        r : float
            Real part of quaternion r + (ix)i + (iy)j + (iz)k
        angle: float
            Angle (in radians) by which to rotate counterclockwise around passed axis
        axis: rebound.Vec3d, list, numpy array, or character 
            3D vector specifying the axis of rotation. The characters "x", "y", "z" are
            shorthand for [1,0,0], [0,1,0], [0,0,1].
        """
        cart = [ix, iy, iz, r]
        angle_axis = [angle, axis]
        if cart.count(None) == len(cart) and angle_axis.count(None) == len(angle_axis):
            super(Rotation, self).__init__(0.0,0.0,0.0,1.0) # Identity
            return
        if cart.count(None) != 0 and angle_axis.count(None) == len(angle_axis):
            raise ValueError("You need to specify all four parameters ix, iy, iz, r.")
        if angle_axis.count(None) != 0 and cart.count(None) == len(cart):
            raise ValueError("You need to specify both angle and axis.")
        if cart.count(None) < len(cart) and angle_axis.count(None) < len(angle_axis):
            raise ValueError("Cannot mix parameters ix, iy, iz, r with angle, axis.")
        if cart.count(None) == 0:
            super(Rotation, self).__init__(ix, iy, iz, r)   
        if angle_axis.count(None) == 0:
            clibrebound.reb_rotation_init_angle_axis.restype = Rotation
            q = clibrebound.reb_rotation_init_angle_axis(c_double(angle), Vec3d(axis)._vec3d)
            super(Rotation, self).__init__(q.ix, q.iy, q.iz, q.r)   

    @classmethod
    def from_to(cls, fromv, tov):
        """
        Returns a Rotation object that maps the 3D fromv vector to the 3D tov vector, i.e., Rotation * fromv = tov.
        Specifically, the rotation is done counterclockwise around the fromv cross tov axis.
        
        Arguments
        ---------
        fromv: list-like (e.g. list, numpy array), or character
            Input 3D vector that Rotation will map to vector tov
        tov: list-like (e.g. list, numpy array), or character
            Output 3D vector when Rotation is applied to fromv
        
        Examples
        --------

        >>> fromv = [2,3,-5]                                            # An arbitrary vector
        >>> rot = rebound.Rotation.from_to(fromv=fromv, tov=[0,0,1])    # Get a rotation that maps fromv to the z axis
        >>> print(rot * fromv)                                          # When our rotation acts on fromv, we get back [0,0,1]
        """
        # "from" is a keyword, need to use somethign else: "fromv"
        _from = Vec3d(fromv)
        _to = Vec3d(tov)
        clibrebound.reb_rotation_init_from_to.restype = cls
        q = clibrebound.reb_rotation_init_from_to(_from._vec3d, _to._vec3d)
        return q

    @classmethod
    def orbit(cls, Omega=0.0, inc=0.0, omega=0.0):
        """
        Consider an orbit, which in a reference coordinate system has standard orbital angles Omega, inc, omega.
        This returns a rotation object 'rot' that rotates vectors into this inclined orbital plane. The vector [1,0,0]
        is rotated such that it points towards the pericenter of the orbit. Murray & Dermott Eq. 2.121 (left hand side)
        
        Arguments
        ---------
        Omega: float
            Longitude of ascending node (default 0)
        inc: float
            Inclination (default 0)
        omega: float
            Argument of pericenter (default 0)

        Examples
        --------

        >>> a, e, inc, Omega, omega = 1, 0.1, 0.2, 0.3, 0.4             # Make up arbitrary numbers
        >>> rot = rebound.Rotation.orbit(Omega=Omega, inc=inc, omega=omega) # Initialize a Rotation to that specific orbit 
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1)
        >>> sim.add(a=a, e=e, inc=inc, Omega=Omega, omega=omega)        # 3D orbit
        >>> sim.add(a=a, e=e)                                           # Orbit in the xy plane
        >>> print(sim.particles[1].xyz)
        >>> print(rot * sim.particles[2].xyz)                           # Same location as other particle after applying rot
        """
        clibrebound.reb_rotation_init_orbit.restype = cls
        q = clibrebound.reb_rotation_init_orbit(c_double(Omega), c_double(inc), c_double(omega))
        return q
    
    @classmethod
    def to_new_axes(cls, newz, newx=None):
        """
        Returns a rotation object that rotates vectors into a new coordinate system with the z axis pointing along 
        the vector newz. If newx is passed, the new x axis will point along newx. If not, the new x axis defaults
        to the intersection betweeen the new xy plane (normal to newz) and the original xy plane (specifically along
        the z cross newz direction). The new y axis completes a right-handed coordinate system. 

        If passed, newx should be perpendicular to newz. This function will only take the component of newx that is 
        perpendicular to newz, in order to avoid any rounding issues. 
        
        Arguments
        ---------
        newz: list-like (e.g. list, numpy array)
            3D vector to use as the new z axis
        newx: list-like (e.g. list, numpy array)
            3D vector to use as the new x axis. If not passed, defaults to the intersection between the new xy plane
            (normal to newz) and the original xy plane (along the z cross newz direction).
        
        Examples
        --------

        >>> sim = rebound.Simulation()                  # Initialize a sim with two Jupiters on arbitrary inclined orbits
        >>> sim.add(m=1)                                                 
        >>> sim.add(m=1e-3, a=1, inc=0.3, Omega=4.2)
        >>> sim.add(m=1e-3, a=2, inc=0.1, Omega=0.5)    # Get a rotation to new axes where z points along total ang. momentum L. Didn't specify newx, 
        >>> rot = rebound.Rotation.to_new_axes(newz=sim.angular_momentum()) # so it points along line of nodes between our original xy plane
        >>> sim.rotate(rot)                             # and the new plane perp. to L. Rotate our simulation into this new reference system
        >>> print(sim.angular_momentum())     # Now the total angular momentum points in the z direction as we expect.
        """
        if not newx: # newx not specified, newx will point along z cross newz (line of nodes)
            clibrebound.reb_vec3d_cross.restype = _Vec3d
            newx = Vec3d(clibrebound.reb_vec3d_cross(Vec3d(0,0,1)._vec3d, Vec3d(newz)._vec3d))
            mag = (newx.x**2 + newx.y**2 + newx.z**2)**(0.5)
            if mag < 1e-15: # z and newz point in the same direction, so line of nodes undefined, don't rotate x
                newx.x, newx.y, newx.z = 1,0,0
        clibrebound.reb_rotation_init_to_new_axes.restype = cls
        q = clibrebound.reb_rotation_init_to_new_axes(Vec3d(newz)._vec3d, Vec3d(newx)._vec3d)
        return q

    def inverse(self):
        """
        Returns a Rotation object's inverse rotation.
        """
        clibrebound.reb_rotation_inverse.restype = Rotation 
        q = clibrebound.reb_rotation_inverse(self)
        return q
    
    def __mul__(self, other):
        if isinstance(other, Rotation):
            clibrebound.reb_rotation_mul.restype = Rotation 
            q = clibrebound.reb_rotation_mul(self, other)
            return q
        if isinstance(other, Particle):
            p = other.copy()
            p.rotate(self)
            return p
        if isinstance(other, Simulation):
            s = other.copy()
            s.rotate(self)
            return s
        try:
            vec = Vec3d(other) # make copy if vec3d, try to convert to vec3d if list-like
            clibrebound.reb_vec3d_irotate(byref(vec._vec3d), self) # rotate vector in place
            return vec
        except:
            return NotImplemented

    def __repr__(self):
        return '<{0}.{1} object at {2}, ix={3}, iy={4}, iz={5}, r={6}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.ix, self.iy, self.iz, self.r)
    _fields_ = [("ix", c_double),
                ("iy", c_double),
                ("iz", c_double),
                ("r", c_double)]

class _Vec3d(Structure):
    """
    Internal use only. Not used as Vec3d directly because assigments to numpy arrays don't worl
    """
    _fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double)]

class Vec3d:
    """
    Class for 3D Cartesian vectors. 
    """
    _vec3d = None

    @property
    def __array_interface__(self):
        return {"shape": (3,), "typestr": "<f8", "data": self._vec3d} 

    def __init__(self, *args):
        if len(args) == 1: 
            vec = args[0]
            if isinstance(vec,_Vec3d):
                vec = [vec.x, vec.y, vec.z]
            elif isinstance(vec,str):
                vec = vec.lower()
                if vec != "x" and vec !="y" and vec != "z":
                    raise ValueError("When passing a string to create a Vec3D, it needs to be one of 'x', 'y', or  'z'")
                if vec == "x":
                    vec = [1,0,0]
                elif vec == "y":
                    vec = [0,1,0]
                if vec == "z":
                    vec = [0,0,1]
            else:
                vec = [float(vec[0]), float(vec[1]), float(vec[2])]
        elif len(args) >= 3:
            vec = [float(args[0]), float(args[1]), float(args[2])]
        self._vec3d =_Vec3d(vec[0],vec[1],vec[2])

    def __mul__(self, other):
        try:
            return Vec3d([self.x*other, self.y*other, self.z*other]) 
        except:
            return NotImplemented

    def __truediv__(self, other):
        if other==0.:
            raise ZeroDivisionError
        try:
            return Vec3d([self.x/other, self.y/other, self.z/other])
        except:
            return NotImplemented

    def __add__(self, other):
        try:
            o = Vec3d(other)
            return Vec3d([self[0]+other[0], self[1]+other[1], self[2]+other[2]]) 
        except:
            return NotImplemented
    
    def __sub__(self, other):
        try:
            o = Vec3d(other)
            return Vec3d([self[0]-other[0], self[1]-other[1], self[2]-other[2]]) 
        except:
            return NotImplemented
    
    def rotate(self, q):
        if not isinstance(q, Rotation):
            raise NotImplementedError
        clibrebound.reb_vec3d_irotate(byref(_vec3d), q)
        return self

    def normalize(self):
        clibrebound.reb_vec3d_normalize.restype = _Vec3d
        r = clibrebound.reb_vec3d_normalize(self._vec3d)
        self._vec3d = r._vec3d
        return self

    def __getitem__(self, key):
        if not isinstance(key, int):
            raise IndexError("Index must be an integer.")
        if key < 0 or key >= 3:
            raise IndexError("Vec3d has exactly three elements and can therefore not access the item with index "+str(key)+".")
        if key == 0:
            return self._vec3d.x
        if key == 1:
            return self._vec3d.y
        if key == 2:
            return self._vec3d.z

    def __setitem__(self, key, value):
        if not isinstance(key, int):
            raise IndexError("Index must be an integer.")
        if key < 0 or key >= 3:
            raise IndexError("Vec3d has exactly three elements and can therefore not access the item with index "+str(key)+".")
        if key == 0:
            self._vec3d.x = c_double(value)
        if key == 1:
            self._vec3d.y = c_double(value)
        if key == 2:
            self._vec3d.z = c_double(value)

    @property
    def x(self):
         return self._vec3d.x
    @x.setter
    def x(self, v):
         self._vec3d.x = v
    @property
    def y(self):
         return self._vec3d.y
    @y.setter
    def y(self, v):
         self._vec3d.y = v
    @property
    def z(self):
         return self._vec3d.z
    @z.setter
    def z(self, v):
         self._vec3d.z = v

    def __repr__(self):
        return '<{0}.{1} object at {2}, [{3}, {4}, {5}]>'.format(self.__module__, type(self).__name__, hex(id(self)), self._vec3d.x, self._vec3d.y, self._vec3d.z)
    

class reb_dp7(Structure):
    _fields_ = [("p0", POINTER(c_double)),
                ("p1", POINTER(c_double)),
                ("p2", POINTER(c_double)),
                ("p3", POINTER(c_double)),
                ("p4", POINTER(c_double)),
                ("p5", POINTER(c_double)),
                ("p6", POINTER(c_double))]

class reb_ghostbox(Structure):
    _fields_ = [("shiftx", c_double),
                ("shifty", c_double),
                ("shiftz", c_double),
                ("shiftvx", c_double),
                ("shiftvy", c_double),
                ("shiftvz", c_double)]

class reb_collision(Structure):
    _fields_ = [("p1", c_int),
                ("p2", c_int),
                ("gb", reb_ghostbox),
                ("ri", c_int)]
    
    def __repr__(self):
        return '<{0}.{1} object at {2}, p1={3}, p2={4}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.p1, self.p2)
    

class reb_simulation_integrator_sei(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_sei.
    It controls the behaviour of the symplectic SEI integrator for shearing
    sheet calculations. It is described in Rein and Tremaine (2011).
    
    This struct should be accessed via the simulation class only. Here is an 
    example:

    >>> sim = rebound.Simulation()
    >>> sim.ri_sei.OMEGA =  1.58
    
    :ivar float OMEGA:          
        The epicyclic frequency OMEGA. For simulations making use of shearing 
        sheet boundary conditions, REBOUND needs to know the epicyclic frequency. 
        By default, OMEGA is 1. For more details read Rein and Tremaine 2011.
    :ivar float OMEGAZ:          
        The z component of the epicyclic frequency OMEGA. By default, it is assuming
        OMEGAZ is the same as OMEGA.
    """
    _fields_ = [("OMEGA", c_double),
                ("OMEGAZ", c_double),
                ("_lastdt", c_double),
                ("_sindt", c_double),
                ("_tandt", c_double),
                ("_sindtz", c_double),
                ("_tandtz", c_double)]

class reb_simulation_integrator_ias15(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_ias15.
    It controls the behaviour of the SEI integrator. See Rein & Spiegel (2015)
    for more information.

    Example usage:

    >>> sim = rebound.Simulation()
    >>> sim.integrator = "ias15"
    >>> sim.ri_ias15.epsilon = 0.

    :ivar float epsilon:          
        Controls the precision of the integrator. Set to 0 for fixed timesteps.
    
    :ivar float min_dt:          
        IAS15 is an adaptive method. This sets the minimum timestep.
    
    :ivar float adaptive_mode:
        Determines how the adaptive timestep is chosen. 
        This replaces the previous epsilon_global variable.
        TODO: list options.
    
    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, epsilon={3}, min_dt={4}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.epsilon, self.min_dt)
    
    _fields_ = [("epsilon", c_double),
                ("min_dt", c_double),
                ("adaptive_mode", c_uint),
                ("_iterations_max_exceeded", c_ulong),
                ("_allocated_N", c_uint),
                ("_at", POINTER(c_double)),
                ("_x0", POINTER(c_double)),
                ("_v0", POINTER(c_double)),
                ("_a0", POINTER(c_double)),
                ("_csx", POINTER(c_double)),
                ("_csv", POINTER(c_double)),
                ("_csa0", POINTER(c_double)),
                ("_g", reb_dp7),
                ("_b", reb_dp7),
                ("_csb", reb_dp7),
                ("_e", reb_dp7),
                ("_br", reb_dp7),
                ("_er", reb_dp7),
                ("_map", POINTER(c_int)),
                ("_map_allocated_n", c_uint),
                ]

class reb_simulation_integrator_saba(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_saba.
    It controls the behaviour of the SABA integrator family.
    See Rein, Tamayo, and Brown (2019) for more details.

    :ivar str type:      
        Set the type of SABA integrator manually. The type can also be set by setting
        the integrator field in the REBOUND simulation. 

    :ivar int safe_mode:      
        This variable acts the same as for WHFast.
        If safe_mode is 1 (default) particles can be modified between
        timesteps and particle velocities and positions are always synchronised.
        If you set safe_mode to 0, the speed and accuracy of the integrator will improve.
        However, make sure you are aware of the consequences. Read the iPython tutorial
        on advanced WHFast usage to learn more.
   
    Example usage:
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "SABA(10,6,4)"
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "SABA"
    >>> sim.ri_saba.type = "(10,6,4)"
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "SABACL4"
    >>> sim.ri_saba.safe_mode = 0

    """
    _fields_ = [("_type", c_uint),
                ("safe_mode", c_uint),
                ("is_synchronized", c_uint),
                ("keep_unsynchronized", c_uint),
            ]
    @property
    def type(self):
        """
        Get or set the type of SABA integrator.
        """
        i = self._type
        for name, _i in SABA_TYPES.items():
            if i==_i:
                return name
        return i
    @type.setter
    def type(self, value):
        if isinstance(value, int):
            self._type = value
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "").replace("(", "").replace(")", "")
            if value in SABA_TYPES: 
                self._type = SABA_TYPES[value]
            else:
                raise ValueError("Warning. SABA type not found.")

class reb_simulation_integrator_whfast(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_whfast.
    It controls the behaviour of the symplectic WHFast integrator described 
    in Rein and Tamayo (2015) and in Rein, Tamayo, and Brown (2019).
    
    This struct should be accessed via the simulation class only. Here is an 
    example:

    >>> sim = rebound.Simulation()
    >>> sim.integrator = "whfast"
    >>> sim.ri_whfast.corrector =  11
    >>> sim.ri_whfast.kernel = "lazy"

    
    :ivar int corrector:      
        The order of the symplectic corrector in the WHFast integrator.
        By default, the symplectic correctors are turned off (=0). For high
        accuracy simulation set this value to 11 or 17. For more details read 
        Rein and Tamayo (2015).
    :ivar int corrector2:      
        Second correctors (C2 of Wisdom et al 1996).
        By default, the second symplectic correctors are turned off (=0). 
        Set to 1 to turn them on.
    :ivar int/string kernel:      
        Kernel option. Set to 0 for the default WH kernel (standard kick step).
        Other options are "modifiedkick" (1), "composition" (2), "lazy" (3).
    :ivar int recalculate_coordinates_this_timestep:
        Sets a flag that tells WHFast that the particles have changed.
        Setting this flag to 1 (default 0) triggers the WHFast integrator
        to recalculate Jacobi/heliocenctric coordinates. This is needed 
        if the user changes the particle position, velocity or mass 
        in-between timesteps.  After every timestep the flag is set back 
        to 0, so if you continuously update the particles manually, 
        you need to set this flag to 1 after every timestep.
    :ivar string coordinates:
        Sets the internal coordinate system that WHFast is using. By default
        it uses ``'jacobi'`` (=0) coordinates. Other options are 
        ``'democraticheliocentric'`` (=1) and ``'whds'`` (=2). See Hernandez 
        and Dehnen (2017) for more information.
    :ivar int safe_mode:
        If safe_mode is 1 (default) particles can be modified between
        timesteps and particle velocities and positions are always synchronised.
        If you set safe_mode to 0, the speed and accuracy of WHFast improve.
        However, make sure you are aware of the consequences. Read the iPython tutorial
        on advanced WHFast usage to learn more.
    """
    _fields_ = [("corrector", c_uint),
                ("corrector2", c_uint),
                ("_kernel", c_uint),
                ("_coordinates", c_uint),
                ("recalculate_coordinates_this_timestep", c_uint),
                ("safe_mode", c_uint),
                ("keep_unsynchronized", c_uint),
                ("_p_jh", POINTER(Particle)),
                ("_p_temp", POINTER(Particle)),
                ("is_synchronized", c_uint),
                ("_allocated_N", c_uint),
                ("_allocated_Ntmp", c_uint),
                ("_timestep_warning", c_uint),
                ("_recalculate_coordinates_but_not_synchronized_warning", c_uint)]

    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, keep_unsynchonized={4}, is_synchronized={5}, corrector={6}, corrector2={7}, kernel={8}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.keep_unsynchronized, self.is_synchronized, self.corrector, self.corrector2, self.kernel)

    @property
    def coordinates(self):
        """
        Get or set the internal coordinate system.

        Available coordinate systems are:

        - ``'jacobi'`` (default)
        - ``'democraticheliocentric'``
        - ``'whds'``
        """
        i = self._coordinates
        for name, _i in WHFAST_COORDINATES.items():
            if i==_i:
                return name
        return i
    @coordinates.setter
    def coordinates(self, value):
        if isinstance(value, int):
            self._coordinates = c_uint(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in WHFAST_COORDINATES: 
                self._coordinates = WHFAST_COORDINATES[value]
            else:
                raise ValueError("Warning. Coordinate system not found.")
    @property
    def kernel(self):
        """
        Get or set the WHFast Kernel.

        Available kernels are:

        - ``'default'`` (standard WH kernel, kick)
        - ``'modifiedkick'`` (modified kick for newtonian gravity)
        - ``'composition'`` (Wisdom's composition method)
        - ``'lazy'`` (Lazy implementer's method)
        """
        i = self._kernel
        for name, _i in WHFAST_KERNELS.items():
            if i==_i:
                return name
        return i
    @kernel.setter
    def kernel(self, value):
        if isinstance(value, int):
            self._kernel = c_uint(value)
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "")
            if value in WHFAST_KERNELS: 
                self._kernel = WHFAST_KERNELS[value]
            else:
                raise ValueError("Warning. Kernel not found.")


class Orbit(Structure):
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
    _fields_ = [("d", c_double),
                ("v", c_double),
                ("h", c_double),
                ("P", c_double),
                ("n", c_double),
                ("a", c_double),
                ("e", c_double),
                ("inc", c_double),
                ("Omega", c_double),
                ("omega", c_double),
                ("pomega", c_double),
                ("f", c_double),
                ("M", c_double),
                ("l", c_double),
                ("theta", c_double),
                ("T", c_double),
                ("rhill", c_double),
                ("pal_h", c_double),
                ("pal_k", c_double),
                ("pal_ix", c_double),
                ("pal_iy", c_double),
                ("hvec", _Vec3d),
                ("evec", _Vec3d)]

    def __str__(self):
        """
        Returns a string with the semi-major axis and eccentricity of the orbit.
        """
        return "<rebound.Orbit instance, a={0} e={1} inc={2} Omega={3} omega={4} f={5}>".format(str(self.a),str(self.e), str(self.inc), str(self.Omega), str(self.omega), str(self.f))
    
    @property 
    def E(self):
        return M_to_E(self.e, self.M)

class Simulation(Structure):
    """
    This is the REBOUND Simulation Class.
    In encapsulated an entire REBOUND simulation and is an abstraction of the C struct reb_simulation.

    ### Examples
    
    Most simulation parameters can be directly changed with the property syntax:

    >>> sim = rebound.Simulation()
    >>> sim.G = 1.                  # Sets the graviational constant (default 1)
    >>> sim.softening = 1.          # Sets the graviational softening parameter (default 0)
    >>> sim.testparticle_type = 1   # Allows massive particles to feel influence from testparticles (default 0)
    >>> sim.dt = 0.1                # Sets the timestep (will change for adaptive integrators such as IAS15).
    >>> sim.t = 0.                  # Sets the current simulation time (default 0)
    >>> print(sim.N)                # Gets the current number of particles
    >>> print(sim.N_active)         # Gets the current number of active particles
       
    By calling rebound.Simulation() as shown above, you create a new simulation object
    The following example creates a simulation, saves it to a file and then creates
    a copy of the simulation store in the binary file.

    >>> sim = rebound.Simulation()
    >>> sim.add(m=1.)
    >>> sim.add(m=1.e-3,x=1.,vy=1.)
    >>> sim.save("simulation.bin")
    >>> sim_copy = rebound.Simulation("simulation.bin")
    
    Similarly, you can create a simulation, from a simulation archive
    by specifying the snapshot you want to load. 

    >>> sim = rebound.Simulation("archive.bin", 34)
    
    or 
    
    >>> sim = rebound.Simulation(filename="archive.bin", snapshot=34)

    Finally, you can also create a new Simulation by passing a bytes object of a 
    SimulationArchive to Simulation():
    
    >>> sim = rebound.Simulation(open("archive.bin","rb").read())

    """
    def __new__(cls, *args, **kw):
        # Handle arguments
        filename = None
        if len(args)>0:
            # If first argument is of type bytes, then this is unpickling a Simulation
            if isinstance(args[0], bytes):
                l = len(args[0]) 
                buft = c_char * l
                buf = buft.from_buffer_copy(args[0])
                # Note: Not calling SimulationArchive.
                # Doing this manually here because we need to keep the reference to buf.
                # So we can later access the contents of the archive to get the simulation.
                sa = SimulationArchive.__new__(SimulationArchive, None, None)
                w = c_int(0)
                clibrebound.reb_read_simulationarchive_from_buffer_with_messages(byref(sa), byref(buf), c_size_t(l), None, byref(w))
                sim = super(Simulation,cls).__new__(cls)
                clibrebound.reb_init_simulation(byref(sim))
                clibrebound.reb_create_simulation_from_simulationarchive_with_messages(byref(sim),byref(sa),c_int(-1),byref(w))
                for majorerror, value, message in BINARY_WARNINGS:
                    if w.value & value:
                        if majorerror:
                            raise RuntimeError(message)
                        else:  
                            # Just a warning
                            warnings.warn(message, RuntimeWarning)
                return sim
            # Otherwise assume first argument is filename
            filename = args[0]
        if "filename" in kw:
            filename = kw["filename"]
        snapshot = -1
        if len(args)>1:
            snapshot = args[1]
        if "snapshot" in kw:
            snapshot = kw["snapshot"]

        # Create simulation
        if filename is None:
            # Create a new simulation
            sim = super(Simulation,cls).__new__(cls)
            clibrebound.reb_init_simulation(byref(sim))
            return sim
        else:
            # Recreate exisitng simulation 
            sa = SimulationArchive(filename,process_warnings=False)
            sim = super(Simulation,cls).__new__(cls)
            clibrebound.reb_init_simulation(byref(sim))
            w = sa.warnings # warnings will be appended to previous warnings (as to not repeat them) 
            clibrebound.reb_create_simulation_from_simulationarchive_with_messages(byref(sim),byref(sa),c_int(snapshot),byref(w))
            for majorerror, value, message in BINARY_WARNINGS:
                if w.value & value:
                    if majorerror:
                        raise RuntimeError(message)
                    else:  
                        # Just a warning
                        warnings.warn(message, RuntimeWarning)
            return sim

    def __init__(self,filename=None,snapshot=None):
        self.save_messages = 1 # Warnings will be checked within python

    def __repr__(self):
        return '<{0}.{1} object at {2}, N={3}, t={4}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.N, self.t)
    
    @classmethod
    def from_archive(cls, filename,snapshot=-1):
        """ rebound.Simulation.from_archive(filename,snapshot) is deprecated and will be removed in the futute. Use rebound.Simulation(filename,snapshot) instead """
        warnings.warn( "rebound.Simulation.from_archive(filename,snapshot) is deprecated and will be removed in the future. Use rebound.Simulation(filename,snapshot) instead", FutureWarning)
        return cls(filename=filename,snapshot=snapshot)

    @classmethod
    def from_file(cls, filename):
        """ rebound.Simulation.from_file(filename) is deprecated and will be removed in the future. Use rebound.Simulation(filename) instead """
        warnings.warn( "rebound.Simulation.from_file(filename) is deprecated and will be removed in the future. Use rebound.Simulation(filename) instead", FutureWarning)
        return cls(filename=filename)
    
    def copy(self):
        """
        Returns a deep copy of a REBOUND simulation. You need to reset 
        any function pointers on the copy. 
        
        Returns
        ------- 
        A rebound.Simulation object.
        
        """
        w = c_int(0)
        sim = Simulation()
        clibrebound.reb_copy_simulation_with_messages(byref(sim),byref(self),byref(w))
        for majorerror, value, message in BINARY_WARNINGS:
            if w.value & value:
                if majorerror:
                    raise RuntimeError(message)
                else:  
                    # Just a warning
                    warnings.warn(message, RuntimeWarning)
        return sim

    def cite(self):
        """
        Generate citations

        This function generates citations to papers relevant to the current 
        setting of the simulation.
        """

        txt, bib = cite(self)
        # one could check for REBOUNDx here, then append txt and bib accordingly
        print(txt + "\n\n\n" + bib)

    def getWidget(self,**kwargs):
        """
        Wrapper function that returns a new widget attached to this simulation.

        Widgets provide real-time 3D visualizations from within a Jupyter notebook.
        See the Widget class for more details on the possible arguments.
        
        
        Arguments
        ---------
        All arguments passed to this wrapper function will be passed to /Widget class.
        
        Returns
        ------- 
        A rebound.Widget object.
        
        Examples
        --------

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(m=1.e-3,x=1.,vy=1.)
        >>> sim.getWidget()

        """
        from .widget import Widget # ondemand
        from ipywidgets import DOMWidget
        from IPython.display import display, HTML
        if not hasattr(self, '_widgets'):
            self._widgets = []
            def display_heartbeat(simp):
                for w in self._widgets:
                    w.refresh(simp,isauto=1)
            self.visualization = VISUALIZATIONS["webgl"] 
            clibrebound.reb_display_init_data(byref(self))
            self._dhbf = AFF(display_heartbeat)
            self._display_heartbeat = self._dhbf
            display(HTML(Widget.getClientCode())) # HACK! Javascript should go into custom.js
        newWidget = Widget(self,**kwargs)
        self._widgets.append(newWidget)
        newWidget.refresh(isauto=0)
        return newWidget
    
    def refreshWidgets(self):
        """
        This function manually refreshed all widgets attached to this simulation.
        
        You want to call this function if any particle data has been manually changed.
        """
        if hasattr(self, '_widgets'):
            for w in self._widgets:
                w.refresh(isauto=0)
        else:
            raise RuntimeError("No widgets found")


# Simulation Archive tools
    def automateSimulationArchive(self, filename, interval=None, walltime=None, step=None, deletefile=False):
        """
        This function automates taking snapshots during a simulation using the Simulation Archive.
        Instead of using this function, one can also call simulationarchive_snapshot() manually
        to create snapshots.

        
        Arguments
        ---------
        filename : str
            Filename of the binary file.
        interval : float
            Interval between outputs in code units.
        walltime : float
            Interval between outputs in wall time (seconds). 
            Useful when using IAS15 with adaptive timesteps. 
        step : int
            Interval between outputs in number of timesteps. 
            Useful when outputs need to be spaced exactly.
        deletefile: bool
            False (default) appends to archive if it exists.
            True deletes filename and starts a new archive.
        
        Examples
        --------
        The following example creates a simulation, then 
        initializes the Simulation Archive and integrates
        it forward in time. 

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(m=1.e-3,x=1.,vy=1.)
        >>> sim.automateSimulationArchive("sa.bin",interval=1000.)
        >>> sim.integrate(1e8)

        The SimulationArchive can later be read in using the following syntax:

        >>> sa = rebound.SimulationArchive("sa.bin")
        >>> sim = sa[0]   # get the first snapshot in the SA file (initial conditions)
        >>> sim = sa[-1]  # get the last snapshot in the SA file

        """
        modes = sum(1 for i in [interval, walltime,step] if i is not None)
        if modes != 1:
            raise AttributeError("Need to specify either interval, walltime, or step")
        if deletefile and os.path.isfile(filename):
            os.remove(filename)
            
            # reset intervals so that automate functions C set sim->next consistently
            self.simulationarchive_auto_interval=0
            self.simulationarchive_auto_walltime=0
            self.simulationarchive_auto_step=0
        if interval:
            clibrebound.reb_simulationarchive_automate_interval(byref(self), c_char_p(filename.encode("ascii")), c_double(interval))
        if walltime:
            clibrebound.reb_simulationarchive_automate_walltime(byref(self), c_char_p(filename.encode("ascii")), c_double(walltime))
        if step:
            clibrebound.reb_simulationarchive_automate_step(byref(self), c_char_p(filename.encode("ascii")), c_ulonglong(step))
        self.process_messages()

    def simulationarchive_snapshot(self, filename, deletefile=False):
        """
        Take a snapshot and save it to a SimulationArchive file.
        If the file does not exist yet, a new one will be created. 
        If the file does exist, a snapshot will be appended.
        
        Arguments
        ---------
        filename : str
            Filename of the binary file.
        deletefile: bool
            False (default) appends to archive if it exists.
            True deletes filename and starts a new archive.

        """
        if deletefile and os.path.isfile(filename):
            os.remove(filename)
        clibrebound.reb_simulationarchive_snapshot(byref(self), c_char_p(filename.encode("ascii")))
        self.process_messages()

    @property
    def simulationarchive_filename(self):
        """
        Returns the current SimulationArchive filename in use. 
        Do not set manually. Use sim.automateSimulationArchive() instead
        """
        return self._simulationarchive_filename

# Message and memory management functions
    def process_messages(self):
        clibrebound.reb_get_next_message.restype = c_int
        buf = create_string_buffer(c_int.in_dll(clibrebound, "reb_max_messages_length").value)
        while clibrebound.reb_get_next_message(byref(self), buf):
            msg = buf.value.decode("ascii")
            if msg[0]=='w':
                warnings.warn(msg[1:], RuntimeWarning)
            elif msg[0]=='e':
                raise RuntimeError(msg[1:])

# Pickling methods: return SimulationArchive binary
    def __reduce__(self):
        buf = c_char_p()
        size = c_size_t()
        clibrebound.reb_output_binary_to_stream(byref(self), byref(buf), byref(size))
        s = bytes(string_at(buf, size=size.value)) #make copy
        clibrebound.reb_output_free_stream(buf) # free original
        return (Simulation, (s,))

# Other operators

    def __del__(self):
        if self._b_needsfree_ == 1: # to avoid, e.g., sim.particles[1]._sim.contents.G creating a Simulation instance to get G, and then freeing the C simulation when it immediately goes out of scope
            clibrebound.reb_free_pointers(byref(self))

    def __eq__(self, other):
        # This ignores the walltime parameter
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_diff_simulations.restype = c_int
        ret = clibrebound.reb_diff_simulations(byref(self), byref(other),c_int(2))
        return not ret
            
    def diff(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_diff_simulations_char.restype = c_char_p
        output = clibrebound.reb_diff_simulations_char(byref(other), byref(self))
        print(output.decode("utf-8"))

    def __add__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        c = self.copy()
        return c.__iadd__(other)
    
    def __iadd__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_simulation_iadd.restype = c_int
        ret = clibrebound.reb_simulation_iadd(byref(self), byref(other))
        if ret==-1:
            raise RuntimeError("Cannot add simulations. Check that the simulations have the same number of particles")
        return self
    
    def __sub__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        c = self.copy()
        return c.__isub__(other)

    def __isub__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_simulation_isub.restype = c_int
        ret = clibrebound.reb_simulation_isub(byref(self), byref(other))
        if ret==-1:
            raise RuntimeError("Cannot subtract simulations. Check that the simulations have the same number of particles")
        return self
    
    def __mul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        c = self.copy()
        c.multiply(other, other)
        return c
    
    def __imul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        self.multiply(other, other)
        return self

    def __rmul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        c = self.copy()
        c.multiply(other, other)
        return c
    
    def __div__(self, other):
        return self.__truediv__(other)
    
    def __idiv__(self, other):
        return self.__itruediv__(other)

    def __truediv__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        c = self.copy()
        if other==0.:
            raise ZeroDivisionError
        c.multiply(1./other, 1./other)
        return c

    def __itruediv__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        if other==0.:
            raise ZeroDivisionError
        self.multiply(1./other, 1./other)
        return self

    def multiply(self, scalar_pos, scalar_vel):
        try:
            scalar_pos = float(scalar_pos)
            scalar_vel = float(scalar_vel)
        except:
            raise ValueError("Cannot multiply simulation with non-scalars.")
        clibrebound.reb_simulation_imul(byref(self), c_double(scalar_pos), c_double(scalar_vel))
    
    def rotate(self, q):
        if not isinstance(q, Rotation):
            return NotImplemented
        clibrebound.reb_simulation_irotate(byref(self), q)

#ODE functions
    def create_ode(self, length, needs_nbody=True):
        clibrebound.reb_create_ode.restype = POINTER(ODE)
        ode_p = clibrebound.reb_create_ode(byref(self), c_int(length))
        ode_p.contents.needs_nbody = c_int(needs_nbody)
        return ODE.from_address(ctypes.addressof(ode_p.contents))

# Status functions
    def status(self, showParticles=True, showAllFields=True):
        """ 
        Prints a summary of the current status 
        of the simulation.
        """
        from rebound import __version__, __build__
        s= ""
        s += "---------------------------------\n"
        s += "REBOUND version:     \t%s\n" %__version__
        s += "REBOUND built on:    \t%s\n" %__build__
        s += "Number of particles: \t%d\n" %self.N       
        s += "Selected integrator: \t" + self.integrator + "\n"       
        s += "Simulation time:     \t%.16e\n" %self.t
        s += "Current timestep:    \t%f\n" %self.dt
        s += "---------------------------------\n"
        if self.N>0 and showParticles:
            for p in self.particles:
                s += str(p) + "\n"
            s += "---------------------------------\n"
        print(s, end="")
        if showAllFields:
            print("The following fields have non-default values:")
            newsim = Simulation()
            clibrebound.reb_diff_simulations_char.restype = c_char_p
            output = clibrebound.reb_diff_simulations_char(byref(newsim), byref(self))
            print(output.decode("utf-8"))



# Set function pointer for additional forces
    @property
    def additional_forces(self):
        """
        Get or set a function pointer for calculating additional forces in the simulation.

        The argument can be a python function or something that can 
        be cast to a C function of type CFUNCTYPE(None,POINTER(Simulaton)). 
        If the forces are velocity dependent, the flag 
        force_is_velocity_dependent needs to be set to 1. Otherwise,
        the particle structures might contain incorrect velocity 
        values.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @additional_forces.setter
    def additional_forces(self, func):
        self._afp = AFF(func)
        self._additional_forces = self._afp

    @property
    def pre_timestep_modifications(self):
        """
        Get or set a function pointer for pre-timestep modifications.

        The argument can be a python function or something that can be cast to a C function or a
        python function.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @pre_timestep_modifications.setter
    def pre_timestep_modifications(self, func):
        self._pretmp = AFF(func)
        self._pre_timestep_modifications = self._pretmp
    
    @property
    def post_timestep_modifications(self):
        """
        Get or set a function pointer for post-timestep modifications.

        The argument can be a python function or something that can be cast to a C function or a
        python function.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @post_timestep_modifications.setter
    def post_timestep_modifications(self, func):
        self._posttmp = AFF(func)
        self._post_timestep_modifications = self._posttmp
 
    @property
    def heartbeat(self):
        """
        Set a function pointer for a heartbeat function.
        The heartbeat function is called every timestep and can be used
        to monitor long simulations, check for stalled simulations and 
        output debugging information.
     
        The argument can be a python function or something that can be 
        cast to a C function or a python function.

        The function called will receive a pointer to the simulation
        object as its argument.
        
        Examples
        --------

        >>> import rebound
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(m=1e-3, a=1)
        >>> def heartbeat(sim):
        >>>     # sim is a pointer to the simulation object,
        >>>     # thus use contents to access object data.
        >>>     # See ctypes documentation for details.
        >>>     print(sim.contents.t)
        >>> sim.heartbeat = heartbeat
        >>> sim.integrate(1.)

        """
        raise AttributeError("You can only set C function pointers from python (not get).")
    @heartbeat.setter
    def heartbeat(self, func):
        self._hb = AFF(func)
        self._heartbeat = self._hb

    @property 
    def coefficient_of_restitution(self):
        """
        Get or set a function pointer that defined the coefficient of restitution.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @coefficient_of_restitution.setter
    def coefficient_of_restitution(self, func):
        self._corfp = CORFF(func)
        self._coefficient_of_restitution = self._corfp
    
    @property 
    def collision_resolve(self):
        """
        Get or set a function pointer for collision resolving routine.
        
        Possible options for setting:
          1) Function pointer
          2) "merge": two colliding particles will merge) 
          3) "hardsphere": two colliding particles will bounce off using a set coefficient of restitution
        """
        raise AttributeError("You can only set C function pointers from python.")
    @collision_resolve.setter
    def collision_resolve(self, func):
        if func == "merge":
            clibrebound.reb_set_collision_resolve.restype = None
            clibrebound.reb_set_collision_resolve(byref(self), clibrebound.reb_collision_resolve_merge)
        elif func == "hardsphere":
            clibrebound.reb_set_collision_resolve.restype = None
            clibrebound.reb_set_collision_resolve(byref(self), clibrebound.reb_collision_resolve_hardsphere)
        elif func == "halt":
            clibrebound.reb_set_collision_resolve.restype = None
            clibrebound.reb_set_collision_resolve(byref(self), clibrebound.reb_collision_resolve_halt)
        else:
            self._colrfp = COLRFF(func)
            self._collision_resolve = self._colrfp
    
    @property 
    def free_particle_ap(self):
        """
        Get or set a function pointer for freeing the ap pointer whenever sim.remove is called.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @free_particle_ap.setter
    def free_particle_ap(self, func):
        self._fpa = FPA(func)
        self._free_particle_ap = self._fpa

# Setter/getter of parameters and constants
    @property 
    def N_real(self):
        """
        Get the current number of real (i.e. no variational/shadow) particles in the simulation.
        """
        return self.N-self.N_var

    @property
    def integrator(self):
        """
        Get or set the integrator module.

        Available integrators include:

        - ``'IAS15'`` (default)
        - ``'WHFast'``
        - ``'SEI'``
        - ``'LEAPFROG'``
        - ``'JANUS'``
        - ``'MERCURIUS'``
        - ``'WHCKL'`` 
        - ``'WHCKM'`` 
        - ``'WHCKC'`` 
        - ``'SABA4'`` 
        - ``'SABACL4'`` 
        - ``'SABACM4'`` 
        - ``'SABA(10,6,4)'`` 
        - ``'EOS'`` 
        - ``'BS'`` 
        - ``'WHFast512'``
        - ``'none'``
        
        Check the online documentation for a full description of each of the integrators. 
        """
        i = self._integrator
        for name, _i in INTEGRATORS.items():
            if i==_i:
                return name
        return i
    @integrator.setter
    def integrator(self, value):
        if isinstance(value, int):
            self._integrator = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in INTEGRATORS: 
                self._integrator = INTEGRATORS[value]
            # Shortcuts
            elif value=="wh":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 0
                self.ri_whfast.kernel = "default"
            elif value=="whc":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "default"
            elif value=="whckl":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "lazy"
            elif value=="whckm":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "modifiedkick"
            elif value=="whckc":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "composition"
            elif value[0:4]=="saba" and len(value)>4:
                self.integrator = "saba"
                self.ri_saba.type = value[4:]
            else:
                raise ValueError("Integrator not found.")
    
    @property
    def boundary(self):
        """
        Get or set the boundary module.

        Available boundary modules are:

        - ``'none'`` (default)
        - ``'open'``
        - ``'periodic'``
        - ``'shear'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._boundary
        for name, _i in BOUNDARIES.items():
            if i==_i:
                return name
        return i
    @boundary.setter
    def boundary(self, value):
        if isinstance(value, int):
            self._boundary = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in BOUNDARIES: 
                self._boundary = BOUNDARIES[value]
            else:
                raise ValueError("Warning. Boundary condition module not found.")

    @property
    def gravity(self):
        """
        Get or set the gravity module.

        Available gravity modules are:

        - ``'none'`` 
        - ``'basic'`` (default)
        - ``'compensated'``
        - ``'tree'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._gravity
        for name, _i in GRAVITIES.items():
            if i==_i:
                return name
        return i
    @gravity.setter
    def gravity(self, value):
        if isinstance(value, int):
            self._gravity = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in GRAVITIES: 
                self._gravity = GRAVITIES[value]
            else:
                raise ValueError("Warning. Gravity module not found.")

    @property
    def collision(self):
        """
        Get or set the collision module.

        Available collision modules are:

        - ``'none'`` (default)
        - ``'direct'``
        - ``'tree'``
        - ``'mercurius'`` 
        - ``'direct'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._collision
        for name, _i in COLLISIONS.items():
            if i==_i:
                return name
        return i
    @collision.setter
    def collision(self, value):
        if isinstance(value, int):
            self._collision = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in COLLISIONS: 
                self.collision = COLLISIONS[value]
            else:
                raise ValueError("Warning. Collision module not found.")

# Units

    @property
    def units(self):
        """
        Tuple of the units for length, time and mass.  Can be set in any order, and strings are not case-sensitive.  See ipython_examples/Units.ipynb for more information.  You can check the units' exact values and add Additional units in rebound/rebound/units.py.  Units should be set before adding particles to the simulation (will give error otherwise).

        Currently supported Units 
        -------------------------

        Times:
        Hr          : Hours
        Yr          : Julian years
        Jyr         : Julian years
        Sidereal_yr : Sidereal year
        Yr2pi       : Year divided by 2pi, with year defined as orbital period of planet at 1AU around 1Msun star
        Kyr         : Kiloyears (Julian)
        Myr         : Megayears (Julian)
        Gyr         : Gigayears (Julian)

        Lengths:
        M           : Meters
        Cm          : Centimeters
        Km          : Kilometers
        AU          : Astronomical Units

        Masses:
        Kg          : Kilograms
        Msun        : Solar masses
        Mmercury    : Mercury masses
        Mvenus      : Venus masses
        Mearth      : Earth masses
        Mmars       : Mars masses
        Mjupiter    : Jupiter masses
        Msaturn     : Saturn masses
        Muranus     : Neptune masses
        Mpluto      : Pluto masses
        
        Examples
        --------
        
        >>> sim = rebound.Simulation()
        >>> sim.units = ('yr', 'AU', 'Msun')

        """
        
        return {'length':hash_to_unit(self.python_unit_l), 'mass':hash_to_unit(self.python_unit_m), 'time':hash_to_unit(self.python_unit_t)}

    @units.setter
    def units(self, newunits):
        newunits = check_units(newunits)        
        if self.N>0: # some particles are loaded
            raise AttributeError("Error:  You cannot set the units after populating the particles array.  See ipython_examples/Units.ipynb.")
        self.update_units(newunits) 

    def update_units(self, newunits): 
        clibrebound.reb_hash.restype = c_uint32
        self.python_unit_l = clibrebound.reb_hash(c_char_p(newunits[0].encode("ascii")))
        self.python_unit_t = clibrebound.reb_hash(c_char_p(newunits[1].encode("ascii")))
        self.python_unit_m = clibrebound.reb_hash(c_char_p(newunits[2].encode("ascii")))
        self.G = convert_G(newunits)

    def convert_particle_units(self, *args): 
        """
        Will convert the units for the simulation (i.e. convert G) as well as the particles' cartesian elements.
        Must have set sim.units ahead of calling this function, so REBOUND knows what units to convert from.

        Parameters
        ----------
        3 strings corresponding to units of time, length and mass.  Can be in any order and aren't case-sensitive.        You can add new units to rebound/rebound/units.py
        """
        if self.python_unit_l == 0 or self.python_unit_m == 0 or self.python_unit_t == 0:
            raise AttributeError("Must set sim.units before calling convert_particle_units in order to know what units to convert from.")
        new_l, new_t, new_m = check_units(args)
        for p in self.particles:
            units_convert_particle(p, hash_to_unit(self.python_unit_l), hash_to_unit(self.python_unit_t), hash_to_unit(self.python_unit_m), new_l, new_t, new_m)
        self.update_units((new_l, new_t, new_m))

# Variational Equations
    def add_variation(self,order=1,first_order=None, first_order_2=None, testparticle=-1):
        """ 
        This function adds a set of variational particles to the simulation. 

        If there are N real particles in the simulation, this functions adds N additional variational 
        particles. To see how many particles (real and variational) are in a simulation, use ``'sim.N'``. 
        To see how many variational particles are in a simulation use ``'sim.N_var'``.

        Currently Leapfrog, WHFast and IAS15 support first order variational equations. IAS15 also
        supports second order variational equations.

        Parameters
        ----------
        order : integer, optional
            By default, the function adds a set of first order variational particles to the simulation. Set this flag to 2 for second order.
        first_order : Variation, optional
            Second order variational equations depend on their corresponding first order variational equations. 
            This parameter expects the Variation object corresponding to the first order variational equations. 
        first_order_2 : Variation, optional
            Same as first_order. But allows to set two different indicies to calculated off-diagonal elements. 
            If omitted, then first_order will be used for both first order equations.
        testparticle : int, optional
            If set to a value >= 0, then only one variational particle will be added and be treated as a test particle.
            

        Returns
        -------
        Returns Variation object (a copy--you can only modify it through its particles property or vary method). 
        """
        cur_var_config_N = self.var_config_N
        if order==1:
            index = clibrebound.reb_add_var_1st_order(byref(self),c_int(testparticle))
        elif order==2:
            if first_order is None:
                raise AttributeError("Please specify corresponding first order variational equations when initializing second order variational equations.")
            if first_order_2 is None:
                first_order_2 = first_order
            index = clibrebound.reb_add_var_2nd_order(byref(self),c_int(testparticle),c_int(first_order.index),c_int(first_order_2.index))
        else:
            raise AttributeError("Only variational equations of first and second order are supported.")

        # Need a copy because location of original might shift if more variations added
        s = Variation.from_buffer_copy(self.var_config[cur_var_config_N])

        return s
        
# MEGNO
    def init_megno(self, seed=None):
        """
        This function initialises the chaos indicator MEGNO particles and enables their integration.

        MEGNO is short for Mean Exponential Growth of Nearby orbits. It can be used to test
        if a system is chaotic or not. In the backend, the integrator is integrating an additional set
        of particles using the variational equation. Note that variational equations are better 
        suited for this than shadow particles. MEGNO is currently only supported in the IAS15 
        and WHFast integrators.

        This function also needs to be called if you are interested in the Lyapunov exponent as it is
        calculate with the help of MEGNO. See Rein and Tamayo 2015 for details on the implementation.

        For more information on MEGNO see e.g. https://dx.doi.org/10.1051/0004-6361:20011189
        """
        if seed is None:
            clibrebound.reb_tools_megno_init(byref(self))
        else:
            clibrebound.reb_tools_megno_init_seed(byref(self), c_uint(seed))
    
    def calculate_megno(self):
        """
        Return the current MEGNO value.
        Note that you need to call init_megno() before the start of the simulation.
        """
        if self._calculate_megno==0:
            raise RuntimeError("MEGNO cannot be calculated. Make sure to call init_megno() after adding all particles but before integrating the simulation.")

        clibrebound.reb_tools_calculate_megno.restype = c_double
        return clibrebound.reb_tools_calculate_megno(byref(self))
    
    def calculate_lyapunov(self):
        """
        Return the current Lyapunov Characteristic Number (LCN).
        Note that you need to call init_megno() before the start of the simulation.
        Different definitions of the LCN exist.  Here, we're following Eq 24 of 
        Cincotta and Simo (2000): https://aas.aanda.org/articles/aas/abs/2000/20/h1686/h1686.html.
        To get a timescale (the Lyapunov timescale), take the inverse of this quantity.
        """
        if self._calculate_megno==0:
            raise RuntimeError("Lyapunov Characteristic Number cannot be calculated. Make sure to call init_megno() after adding all particles but before integrating the simulation.")

        clibrebound.reb_tools_calculate_lyapunov.restype = c_double
        return clibrebound.reb_tools_calculate_lyapunov(byref(self))
    
# Particle add function, used to be called particle_add() and add_particle() 
    def add(self, particle=None, **kwargs):   
        """
        Adds a particle to REBOUND. Accepts one of the following:

        1) A single Particle structure.
        2) The particle's mass and a set of cartesian coordinates: m,x,y,z,vx,vy,vz.
        3) The primary as a Particle structure, the particle's mass and a set of orbital elements: primary,m,a,anom,e,omega,inv,Omega,MEAN (see :class:`.Orbit` for the definition of orbital elements).
        4) A name of an object (uses NASA Horizons to look up coordinates)
        5) A list of particles or names.
        """
        if particle is not None:
            if isinstance(particle, Particle):
                if (self.gravity == "tree" or self.collision == "tree") and self.root_size <=0.:
                    raise ValueError("The tree code for gravity and/or collision detection has been selected. However, the simulation box has not been configured yet. You cannot add particles until the the simulation box has a finite size.")

                clibrebound.reb_add(byref(self), particle)
            elif isinstance(particle, list):
                for p in particle:
                    self.add(p, **kwargs)
            elif isinstance(particle,str):
                if self.python_unit_l == 0 or self.python_unit_m == 0 or self.python_unit_t == 0:
                    self.units = ('AU', 'yr2pi', 'Msun')
                    self.G = 1.0
                builtindatasets = ["solar system", "outer solar system"]
                if particle.lower() == "solar system":          # built in test dataset
                    data.add_solar_system(self)
                elif particle.lower() == "outer solar system":  # built in test dataset
                    data.add_outer_solar_system(self)
                else:
                    if "frame" not in kwargs:
                        if hasattr(self, 'default_plane'):
                            kwargs["plane"] = self.default_plane # allow ASSIST to set default plane
                    self.add(horizons.getParticle(particle, **kwargs), hash=particle)
                    units_convert_particle(self.particles[-1], 'km', 's', 'kg', hash_to_unit(self.python_unit_l), hash_to_unit(self.python_unit_t), hash_to_unit(self.python_unit_m))
            else: 
                raise ValueError("Argument passed to add() not supported.")
        else: 
            self.add(Particle(simulation=self, **kwargs))
        if hasattr(self, '_widgets'):
            self._display_heartbeat(pointer(self))

# Particle getter functions
    @property
    def particles(self):
        """
        Returns a Particles object that allows users to access particles like a dictionary using indices, hashes, or strings. 

        The Particles object uses pointers and thus the contents update 
        as the simulation progresses. Note that the pointers could change,
        for example when a particle is added or removed from the simulation. 
        """
        particles = Particles(self)
        return particles

    @particles.deleter
    def particles(self):
        """
        Remove all particles from the simulation
        """
        clibrebound.reb_remove_all(byref(self))
        self.process_messages()

    def remove(self, index=None, hash=None, keepSorted=True):
        """ 
        Removes a particle from the simulation.

        Parameters
        ----------
        index : int, optional
            Specify particle to remove by index.
        hash : c_uint32 or string, optional
            Specifiy particle to remove by hash (if a string is passed, the corresponding hash is calculated).
        keepSorted : bool, optional
            By default, remove preserves the order of particles in the particles array. 
            Might set it to zero in cases with many particles and many removals to speed things up.
        """
        if index is not None:
            clibrebound.reb_remove(byref(self), index, keepSorted)
        if hash is not None:
            hash_types = c_uint32, c_uint, c_ulong
            PY3 = sys.version_info[0] == 3
            if PY3:
                string_types = str,
                int_types = int,
            else:
                string_types = basestring,
                int_types = int, long
            if isinstance(hash, string_types):
                clibrebound.reb_remove_by_hash(byref(self), rebhash(hash), keepSorted)
            elif isinstance(hash, int_types):
                clibrebound.reb_remove_by_hash(byref(self), c_uint32(hash), keepSorted)
            elif isinstance(hash, hash_types):
                clibrebound.reb_remove_by_hash(byref(self), hash, keepSorted)
        if hasattr(self, '_widgets'):
            self._display_heartbeat(pointer(self))

        self.process_messages()

    def particles_ascii(self, prec=8):
        """
        Returns an ASCII string with all particles' masses, radii, positions and velocities.

        Parameters
        ----------
        prec : int, optional
            Number of digits after decimal point. Default 8.
        """
        s = ""
        for p in self.particles:
            s += (("%%.%de "%prec) * 8)%(p.m, p.r, p.x, p.y, p.z, p.vx, p.vy, p.vz) + "\n"
        if len(s):
            s = s[:-1]
        return s
    
    def add_particles_ascii(self, s):
        """
        Adds particles from an ASCII string. 

        Parameters
        ----------
        s : string
            One particle per line. Each line should include particle's mass, radius, position and velocity.
        """
        for l in s.split("\n"):
            r = l.split()
            if len(r):
                try:
                    r = [float(x) for x in r]
                    p = Particle(simulation=self, m=r[0], r=r[1], x=r[2], y=r[3], z=r[4], vx=r[5], vy=r[6], vz=r[7])
                    self.add(p)
                except:
                    raise AttributeError("Each line requires 8 floats corresponding to mass, radius, position (x,y,z) and velocity (x,y,z).")

# Orbit calculation
    def calculate_orbits(self, primary=None, jacobi_masses=False):
        """ 
        Calculate orbital parameters for all particles in the simulation.
        By default, this functions returns the orbits in Jacobi coordinates. 

        If MEGNO is enabled, variational particles will be ignored.

        Parameters
        ----------

        primary     : rebound.Particle, optional
            Set the primary against which to reference the osculating orbit. Default (use Jacobi center of mass).
            For heliocentric coordinates, pass the central object to this parameter. 
        jacobi_masses: bool
            Whether to use jacobi primary mass in orbit calculation. (Default: False)

        Returns
        -------
        Returns an array of Orbits of length N-1.
        """
        orbits = []
       
        if primary is None:
            jacobi = True
            primary = self.particles[0]
            clibrebound.reb_get_com_of_pair.restype = Particle
        else:
            jacobi = False

        for p in self.particles[1:self.N_real]:
            if jacobi_masses is True:
                interior_mass = primary.m
                # orbit conversion uses mu=G*(p.m+primary.m) so set prim.m=Mjac-m so mu=G*Mjac
                primary.m = self.particles[0].m*(p.m + interior_mass)/interior_mass - p.m
                orbits.append(p.calculate_orbit(primary=primary))
                primary.m = interior_mass # back to total mass of interior bodies to update com
            else:
                orbits.append(p.calculate_orbit(primary=primary))
            if jacobi is True: # update com to include current particle for next iteration
                primary = clibrebound.reb_get_com_of_pair(primary, p)

        return orbits

# COM calculation 
    def calculate_com(self, first=0, last=None):
        """
        Returns the center of momentum for all particles in the simulation.

        Parameters
        ----------
        first: int, optional
            If ``first`` is specified, only calculate the center of momentum starting
            from index=``first``.
        last : int or None, optional
            If ``last`` is specified only calculate the center of momentum up to 
            (but excluding) index=``last``.  Same behavior as Python's range function.

        Examples
        --------
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1, x=-20)
        >>> sim.add(m=1, x=-10)
        >>> sim.add(m=1, x=0)
        >>> sim.add(m=1, x=10)
        >>> sim.add(m=1, x=20)
        >>> com = sim.calculate_com()
        >>> com.x
        0.0 
        >>> com = sim.calculate_com(first=2,last=4) # Considers indices 2,3
        >>> com.x
        5.0

        """
        if last is None:
            last = self.N_real
        clibrebound.reb_get_com_range.restype = Particle
        return clibrebound.reb_get_com_range(byref(self), c_int(first), c_int(last))

# Tools
    def serialize_particle_data(self,**kwargs):
        """
        Fast way to access serialized particle data via numpy arrays.

        This function can directly set the values of numpy arrays to
        current particle data. This is significantly faster than accessing
        particle data via `sim.particles` as all the copying is done 
        on the C side. 
        No memory is allocated by this function.
        It expects correctly sized numpy arrays as arguments. The argument
        name indicates what kind of particle data is written to the array. 
        
        Possible argument names are "hash", "m", "r", "xyz", "vxvyvz", and 
        "xyzvxvyvz". The datatype for the "hash" array needs to be uint32. 
        The other arrays expect a datatype of float64. The lengths of 
        "hash", "m", "r" arrays need to be at least sim.N. The lengths of 
        xyz and vxvyvz need to be at least 3*sim.N. The length of
        "xyzvxvyvz" arrays need to be 6*sim.N. Exceptions are raised 
        otherwise.

        Note that this routine is only intended for special use cases
        where speed is an issue. For normal use, it is recommended to
        access particle data via the `sim.particles` array. Be aware of
        potential issues that arise by directly accesing the memory of
        numpy arrays (see numpy documentation for more details).

        Examples
        --------
        This sets an array to the xyz positions of all particles:

        >>> import numpy as np
        >>> a = np.zeros((sim.N,3),dtype="float64")
        >>> sim.serialize_particle_data(xyz=a)
        >>> print(a)

        To get all current radii of particles:

        >>> a = np.zeros(sim.N,dtype="float64")
        >>> sim.serialize_particle_data(r=a)
        >>> print(a)
        
        To get all current radii and hashes of particles:

        >>> a = np.zeros(sim.N,dtype="float64")
        >>> b = np.zeros(sim.N,dtype="uint32")
        >>> sim.serialize_particle_data(r=a,hash=b)
        >>> print(a,b)

        """
        N = self.N
        possible_keys = ["hash","m","r","xyz","vxvyvz","xyzvxvyvz"]
        d = {x:None for x in possible_keys}
        for k,v in kwargs.items():
            if k in d:
                if k == "hash":
                    if v.dtype!= "uint32":
                        raise AttributeError("Expected 'uint32' data type for '%s' array."%k)
                    if v.size<N:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
                else:
                    if v.dtype!= "float64":
                        raise AttributeError("Expected 'float64' data type for %s array."%k)
                    if k in ["xyz", "vxvyvz"]:
                        minsize = 3*N
                    elif k in ["xyzvxvyvz"]:
                        minsize = 6*N
                    else:
                        minsize = N
                    if v.size<minsize:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            else:
                raise AttributeError("Only '%s' are currently supported attributes for serialization." % "', '".join(d.keys()))

        clibrebound.reb_serialize_particle_data(byref(self), d["hash"], d["m"], d["r"], d["xyz"], d["vxvyvz"], d["xyzvxvyvz"])
    
    def set_serialized_particle_data(self,**kwargs):
        """
        Fast way to set serialized particle data via numpy arrays.
        This is the inverse of Simulation.serialize_particle_data()
        and uses the same syntax
        """

        N = self.N
        possible_keys = ["hash","m","r","xyz","vxvyvz","xyzvxvyvz"]
        d = {x:None for x in possible_keys}
        for k,v in kwargs.items():
            if k in d:
                if k == "hash":
                    if v.dtype!= "uint32":
                        raise AttributeError("Expected 'uint32' data type for '%s' array."%k)
                    if v.size<N:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
                else:
                    if v.dtype!= "float64":
                        raise AttributeError("Expected 'float64' data type for %s array."%k)
                    if k in ["xyz", "vxvyvz"]:
                        minsize = 3*N
                    elif k in ["xyzvxvyvz"]:
                        minsize = 6*N
                    else:
                        minsize = N
                    if v.size<minsize:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            else:
                raise AttributeError("Only '%s' are currently supported attributes for serialization." % "', '".join(d.keys()))

        clibrebound.reb_set_serialized_particle_data(byref(self), d["hash"], d["m"], d["r"], d["xyz"], d["vxvyvz"], d["xyzvxvyvz"])

    def move_to_hel(self):
        """
        This function moves all particles in the simulation to the heliocentric frame.
        Note that the simulation will not stay in the heliocentric frame during integrations
        as the heliocentric frame is not an inertial frame.
        """
        clibrebound.reb_move_to_hel(byref(self))
    
    def move_to_com(self):
        """
        This function moves all particles in the simulation to a center of momentum frame.
        In that frame, the center of mass is at the origin and does not move.
        It makes sense to call this function at the beginning of the integration, especially 
        for the high accuracy integrators IAS15 and WHFast.
        """
        clibrebound.reb_move_to_com(byref(self))

    def calculate_energy(self):
        """
        Returns the sum of potential and kinetic energy of all particles in the simulation.
        """
        warnings.warn( "sim.calculate_energy() is deprecated and will be removed in the future. Use sim.energy() instead", FutureWarning)
        clibrebound.reb_tools_energy.restype = c_double
        return clibrebound.reb_tools_energy(byref(self))
    
    def energy(self):
        """
        Returns the sum of potential and kinetic energy of all particles in the simulation.
        """
        clibrebound.reb_tools_energy.restype = c_double
        return clibrebound.reb_tools_energy(byref(self))
   
    def calculate_angular_momentum(self):
        """
        Returns a list of the three (x,y,z) components of the total angular momentum of all particles in the simulation.
        """
        warnings.warn( "sim.calculate_angular_momentum() is deprecated and will be removed in the future. Use sim.angular_momentum() instead", FutureWarning)
        clibrebound.reb_tools_angular_momentum.restype = _Vec3d
        return Vec3d(clibrebound.reb_tools_angular_momentum(byref(self)))
    
    def angular_momentum(self):
        """
        Returns a list of the three (x,y,z) components of the total angular momentum of all particles in the simulation.
        """
        clibrebound.reb_tools_angular_momentum.restype = _Vec3d
        return Vec3d(clibrebound.reb_tools_angular_momentum(byref(self)))

    def configure_box(self, boxsize, root_nx=1, root_ny=1, root_nz=1):
        """
        Initialize the simulation box.

        This function only needs to be called it boundary conditions other than "none"
        are used. In such a case the boxsize must be known and is set with this function.

        Parameters
        ----------
        boxsize : float
            The size of one root box.
        root_nx, root_ny, root_nz : int, optional
            The number of root boxes in each direction. The total size of the simulation box
            will be ``root_nx * boxsize``, ``root_ny * boxsize`` and ``root_nz * boxsize``.
            By default, there will be exactly one root box in each direction.
        """
        clibrebound.reb_configure_box(byref(self), c_double(boxsize), c_int(root_nx), c_int(root_ny), c_int(root_nz))
        return
   
    def configure_ghostboxes(self, nghostx=0, nghosty=0, nghostz=0):
        """
        Initialize the ghost boxes.

        This function only needs to be called it boundary conditions other than "none" or
        "open" are used. In such a case the number of ghostboxes must be known and is set 
        with this function. 
        
        Parameters
        ----------
        nghostx, nghosty, nghostz : int
            The number of ghost boxes in each direction. All values default to 0 (no ghost boxes).
        """
        clibrebound.nghostx = c_int(nghostx)
        clibrebound.nghosty = c_int(nghosty)
        clibrebound.nghostz = c_int(nghostz)
        return


# Input/Output routines
    def save(self, filename):
        """
        Save the entire REBOUND simulation to a binary file.
        """
        clibrebound.reb_output_binary(byref(self), c_char_p(filename.encode("ascii")))

# Integration
    def step(self):
        """
        Perform exactly one integration step with REBOUND. This function is rarely needed.
        Instead, use integrate().
        """
        clibrebound.reb_step(byref(self))
        self.process_messages()
    
    def steps(self, N_steps):
        """
        Perform exactly N_steps integration steps with REBOUND. This function is rarely needed.
        Instead, use integrate().
        """
        clibrebound.reb_steps(byref(self),c_uint(N_steps))
        self.process_messages()

    def integrate(self, tmax, exact_finish_time=1):
        """
        Main integration function. Call this function when you have setup your simulation and want to integrate it forward (or backward) in time. The function might be called many times to integrate the simulation in steps and create outputs in-between steps.
        
        Parameters
        ----------
        tmax : float
            The final time of your simulation. If the current time is 100, and tmax=200, then after the calling the integrate routine, the time has advanced to t=200. If tmax is larger than or equal to the current time, no integration will be performed.
        exact_finish_time: int, optional
            This argument determines whether REBOUND should try to finish at the exact time (tmax) you give it or if it is allowed to overshoot. Overshooting could happen if one starts at t=0, has a timestep of dt=10 and wants to integrate to tmax=25. With ``exact_finish_time=1``, the integrator will choose the last timestep such that t is exactly 25 after the integration, otherwise t=30. Note that changing the timestep does affect the accuracy of symplectic integrators negatively.
        
        Exceptions
        ----------
        Exceptions are thrown when no more particles are left in the simulation or when a generic integration error occured. 
        If you specified exit_min_distance or exit_max_distance, then additional exceptions might thrown for escaping particles or particles that undergo a clos encounter.
        
        Examples
        -------- 
        The typical usage is as follows. Note the use of ``np.linspace`` to create equally spaced outputs. 
        Using ``np.logspace`` can be used to easily produce logarithmically spaced outputs.

        >>> import numpy as np
        >>> for time in np.linspace(0,100.,10): 
        >>>     sim.integrate(time)
        >>>     perform_output(sim)
        
        """
        self.exact_finish_time = c_int(exact_finish_time)
        ret_value = clibrebound.reb_integrate(byref(self), c_double(tmax))
        if ret_value == 1:
            self.process_messages()
            raise SimulationError("An error occured during the integration.")
        if ret_value == 2:
            if self._odes_N>0:
                raise NoParticles("No particles found. Will exit. Use BS integrator to integrate user-defined ODEs without any particles present.");
            else:
                raise NoParticles("No more particles left in simulation.")
        if ret_value == 3:
            raise Encounter("Two particles had a close encounter (d<exit_min_distance).")
        if ret_value == 4:
            raise Escape("A particle escaped (r>exit_max_distance).")
        if ret_value == 5:
            pass # User caused exit. Do not raise error message
        if ret_value == 6:
            raise KeyboardInterrupt
        if ret_value == 7:
            raise Collision("Two particles collided (d < r1+r2)")
        self.process_messages()

    def stop(self):
        """
        Call this function to stop an integration, for example
        from the heartbeat function.
        """
        clibrebound.reb_stop(byref(self))

    def integrator_reset(self):
        """
        Call this function to reset temporary integrator variables
        """
        clibrebound.reb_integrator_reset(byref(self))

    def integrator_synchronize(self):
        """
        Call this function if safe-mode is disabled and you need to synchronize particle positions and velocities between timesteps.
        """
        clibrebound.reb_integrator_synchronize(byref(self))
    
    def tree_update(self):
        """
        Call this function to update the tree structure manually after removing particles.
        """
        clibrebound.reb_tree_update(byref(self))
    
class Variation(Structure):
    """
    REBOUND Variational Configuration Object.

    This object encapsulated the configuration of one set of variational 
    equations in a REBOUND simulation.  It is an abstraction of the 
    C struct reb_variational_configuration.

    None of the fields in this struct should be changed after it has
    been initialized.

    One rebound simulation can include any number of first and second order 
    variational equations.

    Note that variations are only encoded as particles for convenience.  
    A variational particle's position and velocity should be interpreted as a derivative, i.e. how much that position or velocity varies with respect to the first or second-order variation.  
    See ipython_examples/VariationalEquations.ipynb and Rein and Tamayo (2016) for details.

    Examples
    --------

    >>> sim = rebound.Simulation()          # Create a simulation
    >>> sim.add(m=1.)                       # Add a star
    >>> sim.add(m=1.e-3, a=1.)              #     a planet
    >>> var_config = sim.add_variation()    # Add a set of variational equations. 
    >>> var_config.particles[1].x = 1.      # Initialize the variational particle corresponding to the planet
    >>> sim.integrate(100.)                 # Integrate the simulation
    >>> print(var_config.particles[0].vx)   # Print the velocity of the variational particle corresponding to the star
    """
    _fields_ = [
                ("_sim", POINTER(Simulation)),
                ("order", c_int),
                ("index", c_int),
                ("testparticle", c_int),
                ("index_1st_order_a", c_int),
                ("index_1st_order_b", c_int),
                ("_lrescale", c_double)]

    def vary(self, particle_index, variation, variation2=None, primary=None):
        """
        This function can be used to initialize the variational particles that are 
        part of a Variation.
    
        Note that rather than using this convenience function, one can 
        also directly manipulate the particles' coordinate using the following
        syntax:

        >>> var = sim.add_variation()
        >>> var.particles[0].x = 1.

        The ``vary()`` function is useful for initializing variations corresponding to 
        changes in one of the orbital parameters for a particle on a bound 
        Keplerian orbit.

        The function supports both first and second order variations in the following
        classical orbital parameters:
          a, e, inc, omega, Omega, f
        as well as the Pal (2009) coordinates: 
          a, h, k, ix, iy, lambda
        and in both cases the mass m of the particle. The advantage of the Pal coordinate
        system is that all derivatives are well behaved (infinitely differentiable).
        Classical orbital parameters on the other hand exhibit coordinate singularities, 
        for example when e=0.
        
        The following example initializes the variational particles corresponding to a 
        change in the semi-major axis of the particle with index 1:
        
        >>> var = sim.add_variation()
        >>> var.vary(1,"a")

        Parameters
        ----------
        particle_index : int
            The index of the particle that should be varied. The index starts at 0 and runs through N-1. The first particle added to a simulation receives the index 0, the second 1, and the on.
        variation : string
            This parameter determines which orbital parameter is varied. 
        variation2: string, optional
            This is only used for second order variations which can depend on two varying parameters. If omitted, then it is assumed that the parameter variation is variation2.
        primary: Particle, optional
            By default, variational particles are created in the Heliocentric frame. 
            Set this parameter to use any other particles as a primary (e.g. the center of mass).
        """
        if self.order==2 and variation2 is None:
            variation2 = variation
        if self._sim is not None:
            sim = self._sim.contents
            particles = sim.particles
        else:
            raise RuntimeError("Something went wrong. Cannot seem to find simulation corresponding to variation.")
        if self.testparticle >= 0:
            particles[self.index] = Particle(simulation=sim,particle=particles[particle_index], variation=variation, variation2=variation2, primary=primary)
        else:
            particles[self.index + particle_index] = Particle(simulation=sim,particle=particles[particle_index], variation=variation, variation2=variation2, primary=primary)

    @property
    def particles(self):
        """
        Access the variational particles corresponding to this set of variational equations.

        The function returns a list of particles which are sorted in the same way as those in 
        sim.particles

        The particles are pointers and thus can be modified. 

        If there are N real particles, this function will also return a list of N particles (all of which are 
        variational particles). 
        """
        sim = self._sim.contents
        ps = []
        if self.testparticle>=0:
            N = 1
        else:
            N = sim.N-sim.N_var 
        
        ParticleList = Particle*N
        ps = ParticleList.from_address(ctypes.addressof(sim._particles.contents)+self.index*ctypes.sizeof(Particle))
        return ps
    
    @property
    def lrescale(self):
        """
        Access the lrescale parameter. 
        
        This is a property because sim.add_variation() returns a copy of the struct, so need to find up-to-date reb_variational_configuration struct in simulation.
        """
        sim = self._sim.contents
        for i in range(sim.var_config_N):
            if sim.var_config[i].index == self.index:
                return sim.var_config[i]._lrescale
        raise SimulationError("An error occured while trying to find variational struct in simulation.")
    
    @lrescale.setter
    def lrescale(self, value):
        """
        Set the lrescale parameter. 
        
        This is a property because sim.add_variation() returns a copy of the struct, so need to find up-to-date reb_variational_configuration struct in simulation.
        """
        sim = self._sim.contents
        for i in range(sim.var_config_N):
            if sim.var_config[i].index == self.index:
                self._lrescale = c_double(value)
                sim.var_config[i]._lrescale = c_double(value)
                return

        raise SimulationError("An error occured while trying to find variational struct in simulation.")


        

class reb_particle_int(Structure):
    _fields_ = [
                ("x", c_int64),
                ("y", c_int64),
                ("z", c_int64),
                ("vx", c_int64),
                ("vy", c_int64),
                ("vz", c_int64),
                ]

class reb_simulation_integrator_janus(Structure):
    _fields_ = [
                ("scale_pos",c_double),
                ("scale_vel",c_double),
                ("order", c_uint),
                ("recalculate_integer_coordinates_this_timestep", c_uint),
                ("p_int", POINTER(reb_particle_int)),
                ("_allocated_N",c_uint),
                ]

class reb_simulation_integrator_eos(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_eos.
    It controls the behaviour of the Embedded Operator Splitting methods. See Rein (2019) 
    for more details.
    
    :ivar int,string phi0      
        Sets the Phi_0 operator splitting method
    :ivar int,string phi1     
        Sets the Phi_1 operator splitting method
    :ivar int n     
        Sets the number of substeps taken by Phi_1
    :ivar int safe_mode  
        By default, safe_mode is on (1). Set to 0 (off) to combine
        drift step at the beginning and end of the Phi0 integrator steps.

    Example usage:
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "eos"
    >>> sim.ri_eos.phi0 = "LF8_6_4"
    >>> sim.ri_eos.phi1 = "LF8"
    >>> sim.ri_eos.n = 1
    >>> sim.ri_eos.safe_mode = 0

    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, is_synchronized={4}, n={5}, phi0={6}, phi2={7}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.is_synchronized, self.n, self.phi0, self.phi1)

    @property
    def phi0(self):
        """
        Get or set the type of operator splitting type for phi0.
        """
        i = self._phi0
        for name, _i in EOS_TYPES.items():
            if i==_i:
                return name
        return i
    @phi0.setter
    def phi0(self, value):
        if isinstance(value, int):
            self._phi0 = value
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "").replace("(", "").replace(")", "")
            if value in EOS_TYPES: 
                self._phi0 = EOS_TYPES[value]
            else:
                raise ValueError("Warning. EOS type %s not found."%value)
    @property
    def phi1(self):
        """
        Get or set the type of operator splitting type for phi1.
        """
        i = self._phi1
        for name, _i in EOS_TYPES.items():
            if i==_i:
                return name
        return i
    @phi1.setter
    def phi1(self, value):
        if isinstance(value, int):
            self._phi1 = value
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "").replace("(", "").replace(")", "")
            if value in EOS_TYPES: 
                self._phi1 = EOS_TYPES[value]
            else:
                raise ValueError("Warning. EOS type %s not found."%value)
    _fields_ = [
                ("_phi0",c_uint),
                ("_phi1",c_uint),
                ("n",c_uint),
                ("safe_mode",c_uint),
                ("is_synchronized",c_uint),
                ]

class reb_simulation_integrator_mercurius(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_mercurius.
    It controls the behaviour of the MERCURIUS integrator.  See Rein et al. (2019) 
    for more details.
    
    :ivar float hillfac:      
        Switching radius in units of the hill radius.

    Example usage:
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "mercurius"
    >>> sim.ri_mercurius.hillfac = 3.

    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, is_synchronized={4}, hillfac={5}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.is_synchronized, self.hillfac)

    _fields_ = [("_L", CFUNCTYPE(c_double, POINTER(Simulation), c_double, c_double)),
                ("hillfac", c_double),
                ("recalculate_coordinates_this_timestep", c_uint),
                ("recalculate_dcrit_this_timestep", c_uint),
                ("safe_mode", c_uint),
                ("is_synchronized", c_uint),
                ("mode", c_uint),
                ("_encounterN", c_uint),
                ("_encounterNactive", c_uint),
                ("_tponly_encounter", c_uint),
                ("_allocated_N", c_uint),
                ("_allocated_N_additionalforces", c_uint),
                ("_dcrit_allocated_N", c_uint),
                ("_dcrit", POINTER(c_double)),
                ("_particles_backup", POINTER(Particle)),
                ("_particles_backup_additionalforces", POINTER(Particle)),
                ("_encounter_map", POINTER(c_int)),
                ("_com_pos", _Vec3d),
                ("_com_vel", _Vec3d),
                ]
    @property
    def L(self):
        raise AttributeError("You can only set C function pointers from python.")
    @L.setter
    def L(self, func):
        if func == "mercury":
            self._L = cast(clibrebound.reb_integrator_mercurius_L_mercury,MERCURIUSLF)
        elif func == "C4":
            self._L = cast(clibrebound.reb_integrator_mercurius_L_C4,MERCURIUSLF)
        elif func == "C5":
            self._L = cast(clibrebound.reb_integrator_mercurius_L_C5,MERCURIUSLF)
        elif func == "infinity":
            self._L = cast(clibrebound.reb_integrator_mercurius_L_minfinity,MERCURIUSLF)
        else:
            self._Lfp = MERCURIUSLF(func)
            self._L = self._Lfp

class ODE(Structure):
    @property
    def derivatives(self):
        raise AttributeError("You can only set C function pointers from python.")
    @derivatives.setter
    def derivatives(self, func):
        self._dfp = ODEDER(func)
        func.argtypes = self._dfp.argtypes # I do not understand why this is needed
        self._derivatives = self._dfp
    def update_particles(self):
        clibrebound.reb_integrator_bs_update_particles(self.r, None) 


ODE._fields_ = [
                ("length", c_uint),
                ("allocated_N", c_uint),
                ("needs_nbody", c_int),
                ("y", POINTER(c_double)),
                ("_scale", POINTER(c_double)),
                ("_C", POINTER(c_double)),
                ("_D", POINTER(POINTER(c_double))),
                ("_y1", POINTER(c_double)),
                ("_y0Dot", POINTER(c_double)),
                ("_yDot", POINTER(c_double)),
                ("_yTmp", POINTER(c_double)),
                ("_derivatives", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double), c_double)),
                ("_getscale", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double))),
                ("_pre_timestep", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double))),
                ("_post_timestep", CFUNCTYPE(None,POINTER(ODE), POINTER(c_double))),
                ("r", POINTER(Simulation)),
                ("ref", c_void_p),
            ]               

class reb_simulation_integrator_bs(Structure):
    """
    This class is an abstraction of the C-struct reb_simulation_integrator_bs.
    It controls the behaviour of the Gragg-Bulirsch-Stoer integrator.
    """
    _fields_ = [
                ("_nbody_ode", POINTER(ODE)),
                ("_sequence", POINTER(c_int)),
                ("_costPerStep", POINTER(c_int)),
                ("_costPerTimeUnit", POINTER(c_double)),
                ("_optimalStep", POINTER(c_double)),
                ("_coeff", POINTER(c_double)),
                ("eps_abs", c_double),
                ("eps_rel", c_double),
                ("min_dt", c_double),
                ("max_dt", c_double),
                ("dt_proposed", c_double),
                ("_firstOrLastStep", c_int),
                ("_previousRejected", c_int),
                ("_targetIter", c_int),
                ("_user_ode_needs_nbody", c_int),
            ]               

class timeval(Structure):
    _fields_ = [("tv_sec",c_long),("tv_usec",c_long)]

class reb_display_data(Structure):
    _fields_ = [("r", POINTER(Simulation)),
                ("r_copy", POINTER(Simulation)),
                ("particle_data", c_void_p),
                ("orbit_data", c_void_p),
                ("particles_copy", POINTER(Particle)),
                ("p_jh_copy", POINTER(Particle)),
                ("allocated_N", c_ulong),
                ("allocated_N_whfast", c_ulong),
                ("opengl_enabled", c_int),
                ("scale", c_double),
                ("mouse_x", c_double),
                ("mouse_y", c_double),
                ("retina", c_double),
                # ignoring other data (never used)
                ]

Particle._fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double),
                ("vx", c_double),
                ("vy", c_double),
                ("vz", c_double),
                ("ax", c_double),
                ("ay", c_double),
                ("az", c_double),
                ("m", c_double),
                ("r", c_double),
                ("lastcollision", c_double),
                ("c", c_void_p),
                ("_hash", c_uint32),
                ("ap", c_void_p),
                ("_sim", POINTER(Simulation))]

class reb_simulation_integrator_whfast512(Structure):
    _fields_ = [("is_synchronized", c_uint),
                ("keep_unsynchronized", c_uint),
                ("_allocated_N", c_uint),
                ("gr_potential", c_uint),
                ("recalculate_constants", c_uint),
                ("systems_N", c_uint),
                ("_p_jh", POINTER(Particle)),
                ("_p_jh0", Particle*4)]

# Setting up fields after class definition (because of self-reference)
Simulation._fields_ = [
                ("t", c_double),
                ("G", c_double),
                ("softening", c_double),
                ("dt", c_double),
                ("dt_last_done", c_double),
                ("steps_done", c_ulonglong),
                ("N", c_uint),
                ("N_var", c_int),
                ("var_config_N", c_int),
                ("var_config", POINTER(Variation)),
                ("_var_rescale_warning", c_int),
                ("N_active", c_int),
                ("testparticle_type", c_int),
                ("testparticle_hidewarnings", c_int),
                ("_particle_lookup_table", POINTER(reb_hash_pointer_pair)),
                ("hash_ctr", c_int),
                ("N_lookup", c_int),
                ("allocated_N_lookup", c_int),
                ("allocated_N", c_uint),
                ("_particles", POINTER(Particle)),
                ("gravity_cs", POINTER(_Vec3d)),
                ("gravity_cs_allocated_N", c_int),
                ("_tree_root", c_void_p),
                ("_tree_needs_update", c_int),
                ("opening_angle2", c_double),
                ("_status", c_int),
                ("exact_finish_time", c_int),
                ("force_is_velocity_dependent", c_uint),
                ("gravity_ignore", c_uint),
                ("_output_timing_last", c_double),
                ("_display_clock", c_ulong),
                ("save_messages", c_int),
                ("messages", c_void_p),
                ("exit_max_distance", c_double),
                ("exit_min_distance", c_double),
                ("usleep", c_double),
                ("display_data", POINTER(reb_display_data)),
                ("track_energy_offset", c_int),
                ("energy_offset", c_double),
                ("walltime", c_double),
                ("python_unit_t",c_uint32),
                ("python_unit_l",c_uint32),
                ("python_unit_m",c_uint32),
                ("boxsize", _Vec3d),
                ("boxsize_max", c_double),
                ("root_size", c_double),
                ("root_n", c_int),
                ("root_nx", c_int),
                ("root_ny", c_int),
                ("root_nz", c_int),
                ("nghostx", c_int),
                ("nghosty", c_int),
                ("nghostz", c_int),
                ("collision_resolve_keep_sorted", c_int),
                ("collisions", c_void_p),
                ("collisions_allocated_N", c_int),
                ("minimum_collision_velocity", c_double),
                ("collisions_plog", c_double),
                ("max_radius", c_double*2),
                ("collisions_Nlog", c_long),
                ("_calculate_megno", c_int),
                ("_megno_Ys", c_double),
                ("_megno_Yss", c_double),
                ("_megno_cov_Yt", c_double),
                ("_megno_var_t", c_double),
                ("_megno_mean_t", c_double),
                ("_megno_mean_Y", c_double),
                ("_megno_n", c_long),
                ("rand_seed",c_uint),
                ("simulationarchive_version", c_int),
                ("simulationarchive_size_first", c_long),
                ("simulationarchive_size_snapshot", c_long),
                ("simulationarchive_auto_interval", c_double),
                ("simulationarchive_auto_walltime", c_double),
                ("simulationarchive_auto_step", c_ulonglong),
                ("simulationarchive_next", c_double),
                ("simulationarchive_next_step", c_ulonglong),
                ("_simulationarchive_filename", c_char_p),
                ("_visualization", c_int),
                ("_collision", c_int),
                ("_integrator", c_int),
                ("_boundary", c_int),
                ("_gravity", c_int),
                ("ri_sei", reb_simulation_integrator_sei), 
                ("ri_whfast", reb_simulation_integrator_whfast),
                ("ri_whfast512", reb_simulation_integrator_whfast512),
                ("ri_saba", reb_simulation_integrator_saba),
                ("ri_ias15", reb_simulation_integrator_ias15),
                ("ri_mercurius", reb_simulation_integrator_mercurius),
                ("ri_janus", reb_simulation_integrator_janus),
                ("ri_eos", reb_simulation_integrator_eos),
                ("ri_bs", reb_simulation_integrator_bs),
                ("_odes", POINTER(POINTER(ODE))),
                ("_odes_N", c_int),
                ("_odes_allocated_N", c_int),
                ("_odes_warnings", c_int),
                ("_additional_forces", CFUNCTYPE(None,POINTER(Simulation))),
                ("_pre_timestep_modifications", CFUNCTYPE(None,POINTER(Simulation))),
                ("_post_timestep_modifications", CFUNCTYPE(None,POINTER(Simulation))),
                ("_heartbeat", CFUNCTYPE(None,POINTER(Simulation))),
                ("_display_heartbeat", CFUNCTYPE(None,POINTER(Simulation))),
                ("_coefficient_of_restitution", CFUNCTYPE(c_double,POINTER(Simulation), c_double)),
                ("_collision_resolve", CFUNCTYPE(c_int,POINTER(Simulation), reb_collision)),
                ("_free_particle_ap", CFUNCTYPE(None, POINTER(Particle))),
                ("_extras_cleanup", CFUNCTYPE(None, POINTER(Simulation))),
                ("extras", c_void_p),
                 ]


POINTER_REB_SIM = POINTER(Simulation) 
AFF = CFUNCTYPE(None,POINTER_REB_SIM)
ODEDER = CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double), c_double)
ODESCALE = CFUNCTYPE(None,POINTER(ODE), POINTER(c_double), POINTER(c_double))
CORFF = CFUNCTYPE(c_double,POINTER_REB_SIM, c_double)
COLRFF = CFUNCTYPE(c_int, POINTER_REB_SIM, reb_collision)
MERCURIUSLF = CFUNCTYPE(c_double, POINTER_REB_SIM, c_double, c_double)
FPA = CFUNCTYPE(None, POINTER(Particle))

class Particles(MutableMapping):
    """
    This class allows the user to access particles like a dictionary using the particle's 1) index 2) hash 3) string (which will be converted to hash).
    Allows for negative indices and slicing.
    """
    def __init__(self, sim):
        self.sim = sim

    @property
    def _ps(self):
        ParticleList = Particle*self.sim.N
        pl = ParticleList.from_address(ctypes.addressof(self.sim._particles.contents))
        pl._sim = self.sim # keep reference to sim until ParticleList is deallocated to avoid memory issues
        return pl

    def __getitem__(self, key):
        hash_types = c_uint32, c_uint, c_ulong
        PY3 = sys.version_info[0] == 3
        if PY3:
            string_types = str,
            int_types = int,
        else:
            string_types = basestring,
            int_types = int, long,
       
        if isinstance(key, slice):
            return [self[i] for i in range(*key.indices(len(self)))]

        if isinstance(key, int_types):
            if key < 0: # accept negative indices
                key += self.sim.N
            if key < 0 or key >= self.sim.N:
                raise AttributeError("Index {0} used to access particles out of range.".format(key))
            return self._ps[key]

        else:
            clibrebound.reb_get_particle_by_hash.restype = POINTER(Particle)
            if isinstance(key, string_types):
                key = rebhash(key)
            elif not isinstance(key, hash_types):
                raise AttributeError("Expecting string, integer or ctypes.c_uint32 as argument to sim.particles.  See UniquelyIdentifyingParticlesWithHashes.ipynb ipython_example.")

            ptr = clibrebound.reb_get_particle_by_hash(byref(self.sim), key) 
            
            if ptr:
                p = Particle
                return p.from_address(ctypes.addressof(ptr.contents))
            else:
                raise ParticleNotFound("Particle was not found in the simulation.") 

    def __setitem__(self, key, value):
        if isinstance(value, Particle):
            value._sim = pointer(self.sim)
            p = self[key]
            if p.index == -1:
                raise AttributeError("Can't set particle (particle not found in simulation).")
            else:
                self._ps[p.index] = value

    def __delitem__(self, key):
        pass

    def __iter__(self):
        if self.sim.N>0:
            for p in self._ps:
                yield p

    def __len__(self):
        return self.sim.N

# Import at the end to avoid circular dependence
from . import horizons
from . import data
from .simulationarchive import SimulationArchive
