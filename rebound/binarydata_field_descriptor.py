import ctypes
from . import clibrebound, string_size_max
from .vectors import Vec3dBasic

# Note: ignoring some types that are not supposed to be set by the user
REB_BINARYDATA_DTYPE = {1: ctypes.c_double, 2: ctypes.c_int, 3: ctypes.c_uint, 4: ctypes.c_uint32, 5: ctypes.c_int64,
                        6: ctypes.c_uint64, 7: ctypes.c_size_t, 8: Vec3dBasic, 10: ctypes.c_void_p, 11: ctypes.c_void_p}

class EnumDescriptorListIterator:
    def __init__(self, enum_descriptor_list):
        self._enum_descriptor = enum_descriptor_list
        self._index = 0
    def __next__(self):
        if not self._enum_descriptor_list:
            raise StopIteration
        ed = self._enum_descritor_list[self._index]
        if ed.name == b"":
            raise StopIteration
        self._index += 1 
        return ed
    def __iter__(self):
        return self

class EnumDescriptor(ctypes.Structure):
    _fields_ = [("value", ctypes.c_int), 
                ("name", ctypes.c_char*string_size_max),
                ]
BaseEnumDescriptorList = ctypes.POINTER(EnumDescriptor)
class EnumDescriptorList(BaseEnumDescriptorList):
    def __iter__(self):
        return EnumDescriptorListIterator(self)


class BinarydataFieldDescriptor(ctypes.Structure):
    """ 
    See C structure reb_binarydata_field_descriptor.
    Describes the binarydata field in simulationarchives.
    """

    def __repr__(self):
        return '<{0}.{1} object at {2}, dtype={3}, name=\'{4}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.dtype, self.name.decode("ascii"))
    _fields_ = [("documentation", ctypes.c_char_p),
                ("dtype", ctypes.c_int),
                ("name", ctypes.c_char*string_size_max),
                ("offset", ctypes.c_size_t),
                ("offset_N", ctypes.c_size_t),
                ("element_size", ctypes.c_size_t),
                ("enum_descriptor_list", EnumDescriptorList),
                ]


