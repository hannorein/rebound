import ctypes 

class Rotation(ctypes.Structure):
    """
    This class facilitates rotations of Vec3d objects, and provides various convenience functions
    for commonly used rotations in celestial mechanics.
    """
    def __init__(self, ix=None, iy=None, iz=None, r=None, angle=None, axis=None, fromv=None, tov=None):
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
        angleaxis = [angle, axis]
        fromto = [fromv, tov]
        supplied = [a.count(None)!=len(a) for a in [cart, angleaxis, fromto]]
        if sum(supplied) > 1:
            raise ValueError("Cannot mix parameters.")
        if cart.count(None) == len(cart) and angleaxis.count(None) == len(angleaxis) and fromto.count(None) == len(fromto):
            clibrebound.reb_rotation_identity.restype = Rotation
            q = clibrebound.reb_rotation_identity()
            super(Rotation, self).__init__(q.ix, q.iy, q.iz, q.r)   
            return
        if cart.count(None) != 0 and cart.count(None) != len(cart):
            raise ValueError("You need to specify all four parameters ix, iy, iz, r.")
        if angleaxis.count(None) != 0 and angleaxis.count(None) != len(angleaxis):
            raise ValueError("You need to specify both angle and axis.")
        if fromto.count(None) != 0 and fromto.count(None) != len(fromto):
            raise ValueError("You need to specify both fromv and tov.")
        if cart.count(None) == 0:
            super(Rotation, self).__init__(ix, iy, iz, r)   
            return
        if angleaxis.count(None) == 0:
            clibrebound.reb_rotation_init_angle_axis.restype = Rotation
            q = clibrebound.reb_rotation_init_angle_axis(ctypes.c_double(angle), Vec3d(axis)._vec3d)
            super(Rotation, self).__init__(q.ix, q.iy, q.iz, q.r)   
            return
        if fromto.count(None) == 0:
            q = Rotation.from_to(fromv, tov)
            super(Rotation, self).__init__(q.ix, q.iy, q.iz, q.r)   
            return

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
        q = clibrebound.reb_rotation_init_orbit(ctypes.c_double(Omega), ctypes.c_double(inc), ctypes.c_double(omega))
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
            clibrebound.reb_vec3d_cross.restype = Vec3dBasic
            newx = Vec3d(clibrebound.reb_vec3d_cross(Vec3d(0,0,1)._vec3d, Vec3d(newz)._vec3d))
            mag = (newx.x**2 + newx.y**2 + newx.z**2)**(0.5)
            if mag < 1e-15: # z and newz point in the same direction, so line of nodes undefined, don't rotate x
                newx.x, newx.y, newx.z = 1,0,0
        clibrebound.reb_rotation_init_to_new_axes.restype = cls
        q = clibrebound.reb_rotation_init_to_new_axes(Vec3d(newz)._vec3d, Vec3d(newx)._vec3d)
        return q
    
    def orbital(self):
        """
        Returns a three vector with orbital elements Omega, inc, omega.
        Note: the angles might not always be in the correct quadrant and might be
        inconsistent with REBOUND's standard definition of orbital elements.
        """
        Omega = ctypes.c_double()
        inc = ctypes.c_double()
        omega = ctypes.c_double()
        clibrebound.reb_rotation_to_orbital(self, ctypes.byref(Omega), ctypes.byref(inc), ctypes.byref(omega))
        return [Omega.value, inc.value, omega.value]
    

    def inverse(self):
        """
        Returns a Rotation object's inverse rotation.
        """
        clibrebound.reb_rotation_inverse.restype = Rotation 
        q = clibrebound.reb_rotation_inverse(self)
        return q
    
    def normalize(self):
        """
        Returns a normalized copy of the Rotation object.
        """
        clibrebound.reb_rotation_normalize.restype = Rotation 
        q = clibrebound.reb_rotation_normalize(self)
        return q
    
    def __eq__(self, other):
        if not isinstance(other, Rotation):
            return NotImplemented
        return self.ix == other.ix and self.iy == other.iy and self.iz == other.iz and self.r == other.r
    
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
            clibrebound.reb_vec3d_irotate(ctypes.byref(vec._vec3d), self) # rotate vector in place
            return vec
        except:
            return NotImplemented

    def __repr__(self):
        return '<{0}.{1} object at {2}, ix={3}, iy={4}, iz={5}, r={6}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.ix, self.iy, self.iz, self.r)
    _fields_ = [("ix", ctypes.c_double),
                ("iy", ctypes.c_double),
                ("iz", ctypes.c_double),
                ("r", ctypes.c_double)]

from . import clibrebound
from .simulation import Simulation
from .vectors import Vec3d, Vec3dBasic
from .particle import Particle

