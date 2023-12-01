from ctypes import Structure, c_double

class Vec6d(Structure):
    _fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double),
                ("vx", c_double),
                ("vy", c_double),
                ("vz", c_double)]

class Vec3dBasic(Structure):
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
            if isinstance(vec,Vec3dBasic):
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
        self._vec3d =Vec3dBasic(vec[0],vec[1],vec[2])

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
        clibrebound.reb_vec3d_normalize.restype = Vec3dBasic
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
    
