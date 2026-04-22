import ctypes
from . import clibrebound
from .vectors import Vec3dBasic

# Note: ignoring some types that are not supposed to be set by the user
REB_BINARYDATA_DTYPE = {0: ctypes.c_double, 1: ctypes.c_int, 2: ctypes.c_uint, 3: ctypes.c_uint32, 4: ctypes.c_int64,
                        5: ctypes.c_uint64, 7: Vec3dBasic, 9: ctypes.c_void_p, 10: ctypes.c_void_p, 18: ctypes.c_size_t}

class EnumDescriptor(ctypes.Structure):
    _fields_ = [("value", ctypes.c_int), 
                ("name", ctypes.c_char*256),
                ]

class BinaryFieldDescriptor(ctypes.Structure):
    """ 
    Describes the binary field in simulationarchives 

    Used here for unit tests only. Checking if ids are unique.
    """

    def __repr__(self):
        return '<{0}.{1} object at {2}, type={3}, dtype={4}, name=\'{5}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.type, self.dtype, self.name.decode("ascii"))
    _fields_ = [("type", ctypes.c_uint),
                ("dtype", ctypes.c_int),
                ("name", ctypes.c_char*256),
                ("offset", ctypes.c_size_t),
                ("offset_N", ctypes.c_size_t),
                ("element_size", ctypes.c_size_t),
                ("enum_descriptor_list", ctypes.POINTER(EnumDescriptor)),
                ]

class Integrator(ctypes.Structure):
    _fields_ = [("name", ctypes.c_char_p),
                ("step", ctypes.c_void_p),
                ("synchronize", ctypes.c_void_p),
                ("create", ctypes.c_void_p),
                ("free", ctypes.c_void_p),
                ("did_add_particle", ctypes.c_void_p),
                ("will_remove_particle", ctypes.c_void_p),
                ("field_descriptor_list", ctypes.POINTER(BinaryFieldDescriptor)),
                ("state", ctypes.c_void_p),
                ]
    def __str__(self):
        if self.name:
            return self.name.decode("utf-8")
        else:
            return "<no name provided>"
        return None

    def __eq__(self, other):
        return self.__str__() == other.__str__()

    def __getattr__(self, name):
        i=0
        while True:
            field_descriptor = self.field_descriptor_list[i]
            if field_descriptor.type == 0:
                raise AttributeError("Field '%s' not found in IntegratorState." % name)
            if field_descriptor.name.decode("utf-8") == name:
                if field_descriptor.dtype not in REB_BINARYDATA_DTYPE:
                    raise NotImplemented("Datatype of field '%s' in IntegratorState is currently not supported." % name)
                return REB_BINARYDATA_DTYPE[field_descriptor.dtype].from_address(self.state + field_descriptor.offset).value
            i += 1

    def __setattr__(self, name, value):
        i=0
        while True:
            field_descriptor = self.field_descriptor_list[i]
            if field_descriptor.type == 0:
                raise AttributeError("Field '%s' not found in IntegratorState." % name)
            if field_descriptor.name.decode("utf-8") == name:
                if field_descriptor.dtype not in REB_BINARYDATA_DTYPE:
                    raise NotImplemented("Datatype of field '%s' in IntegratorState is currently not supported." % name)
                pointer_to_field = REB_BINARYDATA_DTYPE[field_descriptor.dtype].from_address(self.state + field_descriptor.offset)
                if isinstance(value, str) and field_descriptor.enum_descriptor_list:
                    j=0
                    while True:
                        enum_descriptor = field_descriptor.enum_descriptor_list[j]
                        if enum_descriptor.name == b"": # reached end of list
                            available_options = []
                            k=0
                            while True:
                                enum_descriptor = field_descriptor.enum_descriptor_list[k]
                                k += 1
                                if enum_descriptor.name == b"": # reached end of list
                                    raise AttributeError("Field '%s' can not be set to '%s' Available options are: %s." % (name, value, ", ".join(available_options)))
                                else:
                                    available_options.append("'"+enum_descriptor.name.decode("utf-8").lower()+"'")
                        if enum_descriptor.name.decode("utf-8") == value.upper():
                            pointer_to_field.value = enum_descriptor.value
                            return
                        j += 1
                pointer_to_field.value = value
                return
            i += 1
    
    def __repr__(self):
        if self.name:
            return '<{0}.{1} object at {2}, name=\'{3}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.name.decode("ascii"))
        else:
            return '<{0}.{1} object at {2}, name=None>'.format(self.__module__, type(self).__name__, hex(id(self)))


def binary_field_descriptor_list():
    fd_pointer_t = ctypes.POINTER(BinaryFieldDescriptor)
    fd_pointer = (BinaryFieldDescriptor*3).in_dll(clibrebound, "reb_binarydata_field_descriptor_list")
    fd_pointer = ctypes.cast(fd_pointer, fd_pointer_t) # not sure why I have to do it this way
    l = []
    i=0
    while True:
        l.append(fd_pointer[i])
        if fd_pointer[i].name == b'end':
            break
        i += 1
    return l


# All available integrators
#sym_address = ctypes.addressof(ctypes.c_void_p.in_dll(clibrebound, "reb_integrators_available"))
#integrators_available = ctypes.cast(sym_address, ctypes.POINTER(ctypes.POINTER(Integrator)))

