import ctypes
import textwrap
from .binarydata_field_descriptor import BinarydataFieldDescriptorList, REB_BINARYDATA_DTYPE, format_doc

class Integrator(ctypes.Structure):
    """ Contains function pointers to integrator routines, documentation, and field data """
    _fields_ = [("documentation", ctypes.c_char_p),
                ("step", ctypes.c_void_p),
                ("synchronize", ctypes.c_void_p),
                ("create", ctypes.c_void_p),
                ("free", ctypes.c_void_p),
                ("did_add_particle", ctypes.c_void_p),
                ("will_remove_particle", ctypes.c_void_p),
                ("field_descriptor_list", BinarydataFieldDescriptorList),
                ]

class IntegratorConfiguration(ctypes.Structure):
    """ Generic Integrator Configuration """
    _fields_ = [("name", ctypes.c_char_p),
                ("callbacks", Integrator),
                ("state", ctypes.c_void_p),
                ]
    @property
    def __doc__(self):
        """ Dynamically generate __doc__ from data in C structs """
        callbacks = self.callbacks
        if callbacks:
            doc = "REBOUND Integrator ("+self.name.decode("utf-8").upper()+")\n\n"
            if callbacks.documentation:
                doc += format_doc(callbacks.documentation) + "\n"
            fdlist = callbacks.field_descriptor_list
            fdlist_doc = []
            for fd in fdlist:
                fd_doc = fd.doc()
                if fd_doc: fdlist_doc.append(fd_doc)
            if len(fdlist_doc):
                doc += "\nAttributes\n"
                doc += "----------\n"
                doc += "\n".join(fdlist_doc)
            return doc
        return "No documentation available."

    def __str__(self):
        """ The string representation of an integrator is its name """
        _name = self.name
        if _name:
            return _name.decode("utf-8")
        else:
            return "<no name provided>"
        return None

    def __eq__(self, value):
        """ Integrators are equal if they have the same name """
        return self.__str__() == value.__str__()

    def __getattr__(self, name):
        """ Get the value of a field in the integrator's state """
        field_descriptor = self.callbacks.field_descriptor_list.field_descriptor_for_name(name)
        value = REB_BINARYDATA_DTYPE[field_descriptor.dtype].from_address(self.state + field_descriptor.offset).value
        return field_descriptor.getInstance(value)

    def __setattr__(self, name, value):
        """ Set a field in the integrator's state to a value """
        field_descriptor = self.callbacks.field_descriptor_list.field_descriptor_for_name(name)
        field_descriptor = self.callbacks.field_descriptor_list.field_descriptor_for_name(name)
        pointer_to_field = REB_BINARYDATA_DTYPE[field_descriptor.dtype].from_address(self.state + field_descriptor.offset)
        if isinstance(value, str) and field_descriptor.enum_descriptor_list:
            for enum_descriptor in field_descriptor.enum_descriptor_list:
                if enum_descriptor.name.decode("utf-8").upper() == value.upper():
                    pointer_to_field.value = enum_descriptor.value
                    return
            raise AttributeError("Field '%s' can not be set to '%s' Available options are: %s." % (name, value, ", ".join([ed.name.decode("utf-8") for ed in field_descriptor.enum_descriptor_list])))
        pointer_to_field.value = value
    
    def __repr__(self):
        """ Create a repr with all integrator fields shown for easy debugging """
        fields = {"name": self.name.decode("utf-8")}
        for fd in self.callbacks.field_descriptor_list:
            if fd.dtype in REB_BINARYDATA_DTYPE and REB_BINARYDATA_DTYPE[fd.dtype] is not None:
                value = REB_BINARYDATA_DTYPE[fd.dtype].from_address(self.state + fd.offset).value
                for enum_descriptor in fd.enum_descriptor_list:
                    if enum_descriptor.value == value:
                        value = enum_descriptor.name.decode("utf-8")
                        break
                if value and REB_BINARYDATA_DTYPE[fd.dtype] == ctypes.c_void_p:
                    value = hex(value) # pointers as hex
                fields[fd.name.decode("utf-8")] = value
            else:
                fields[fd.name.decode("utf-8")] = "<value not shown>"
        values = ", ".join([name+"="+str(value) for name,value in fields.items()])
        return '<{0}.{1} object at {2}, {3}>'.format(self.__module__, type(self).__name__, hex(id(self)), values)

