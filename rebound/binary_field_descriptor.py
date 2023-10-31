import ctypes
from . import clibrebound

class BinaryFieldDescriptor(ctypes.Structure):
    """ 
    Describes the binary field in simulationarchives 

    Used here for unit tests only. Checking if ids are unique.
    """

    def __repr__(self):
        return '<{0}.{1} object at {2}, type={3}, dtype={4}, name=\'{5}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.type, self.dtype, self.name.decode("ascii"))
    _fields_ = [("type", ctypes.c_uint),
                ("dtype", ctypes.c_int),
                ("name", ctypes.c_char*1024),
                ("offset", ctypes.c_size_t),
                ("offset_N", ctypes.c_size_t),
                ("element_size", ctypes.c_size_t),
                ]
def binary_field_descriptor_list():
    fd_pointer_t = ctypes.POINTER(BinaryFieldDescriptor)
    fd_pointer = (BinaryFieldDescriptor*3).in_dll(clibrebound, "reb_binary_field_descriptor_list")
    fd_pointer = ctypes.cast(fd_pointer, fd_pointer_t) # not sure why I have to do it this way
    l = []
    i=0
    while True:
        l.append(fd_pointer[i])
        if fd_pointer[i].name == b'end':
            break
        i += 1
    return l


