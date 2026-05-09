import ctypes
from .binarydata_field_descriptor import BinarydataFieldDescriptor, REB_BINARYDATA_DTYPE

class Integrator(ctypes.Structure):
    _fields_ = [("documentation", ctypes.c_char_p),
                ("step", ctypes.c_void_p),
                ("synchronize", ctypes.c_void_p),
                ("create", ctypes.c_void_p),
                ("free", ctypes.c_void_p),
                ("did_add_particle", ctypes.c_void_p),
                ("will_remove_particle", ctypes.c_void_p),
                ("field_descriptor_list", ctypes.POINTER(BinarydataFieldDescriptor)),
                ]

class IntegratorConfiguration(ctypes.Structure):
    _fields_ = [("name", ctypes.c_char_p),
                ("callbacks", Integrator),
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

    def _field_descriptor_for_name(self, name):
        i=0
        while self.callbacks.field_descriptor_list:
            field_descriptor = self.callbacks.field_descriptor_list[i]
            if field_descriptor.name == b'':
                raise AttributeError("Field '%s' not found in IntegratorState." % name)
            if field_descriptor.name.decode("utf-8") == name:
                if field_descriptor.dtype not in REB_BINARYDATA_DTYPE:
                    raise NotImplemented("Datatype of field '%s' in IntegratorState is currently not supported." % name)
                return field_descriptor
            i += 1

    def __getattr__(self, name):
        field_descriptor = self._field_descriptor_for_name(name)
        return REB_BINARYDATA_DTYPE[field_descriptor.dtype].from_address(self.state + field_descriptor.offset).value

    def __setattr__(self, name, value):
        field_descriptor = self._field_descriptor_for_name(name)
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
    
    def __repr__(self):
        if self.name:
            return '<{0}.{1} object at {2}, name=\'{3}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.name.decode("ascii"))
        else:
            return '<{0}.{1} object at {2}, name=None>'.format(self.__module__, type(self).__name__, hex(id(self)))

