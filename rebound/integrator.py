import ctypes
import textwrap
from .binarydata_field_descriptor import BinarydataFieldDescriptorList, REB_BINARYDATA_DTYPE

class Integrator(ctypes.Structure):
    _fields_ = [("documentation", ctypes.c_char_p),
                ("step", ctypes.c_void_p),
                ("synchronize", ctypes.c_void_p),
                ("create", ctypes.c_void_p),
                ("free", ctypes.c_void_p),
                ("did_add_particle", ctypes.c_void_p),
                ("will_remove_particle", ctypes.c_void_p),
                ("field_descriptor_list", BinarydataFieldDescriptorList),
                ]

import re
from textwrap import TextWrapper

class MarkdownTextWrapper(TextWrapper):
    """A TextWrapper which handles  markdowni reference links."""
    
    LINK_REGEX = re.compile(r"(\[.*?\](?::\s*https?://[^\s]+))")

    def _split(self, text):
        split = re.split(self.LINK_REGEX, text)
        chunks: List[str] = []
        for item in split:
            match = re.match(self.LINK_REGEX, item)
            if match:
                chunks.append(item) # do not break links
            else:
                chunks.extend(super()._split(item)) # handle normally
        return chunks

def format_doc(doc, indent=""):
    doc = doc.decode("utf-8")
    lines = doc.splitlines()
    doc = ""
    wrapper = MarkdownTextWrapper(width=80, initial_indent=indent, subsequent_indent=indent, break_long_words=False)
    for i, line in enumerate(lines):
        if i:
            doc += "\n"
        doc += wrapper.fill(line)
    return doc

class IntegratorConfiguration(ctypes.Structure):
    """ Generic Integrator Configuration """
    _fields_ = [("name", ctypes.c_char_p),
                ("callbacks", Integrator),
                ("state", ctypes.c_void_p),
                ]
    # Automatically generate __doc__ from data in C structs
    @property
    def __doc__(self):
        callbacks = self.callbacks
        if callbacks:
            doc = "REBOUND Integrator ("+self.name.decode("utf-8")+")\n\n"
            if callbacks.documentation:
                doc += format_doc(callbacks.documentation) + "\n"
            fdlist = callbacks.field_descriptor_list
            attributes = 0
            for fd in fdlist:
                if fd.name == b'':
                    break
                if fd.documentation == b'':
                    continue
                if attributes == 0:
                    doc += "\nAttributes\n"
                    doc += "----------\n"
                else:
                    doc += "\n"
                attributes += 1
                doc += fd.name.decode("utf-8") + " : " 
                if fd.dtype not in REB_BINARYDATA_DTYPE:
                    doc += "(unknown datatype)\n"
                else:
                    doc += str(REB_BINARYDATA_DTYPE[fd.dtype].__name__) + "\n"
                doc += format_doc(fd.documentation,indent="    ")
                edl = fd.enum_descriptor_list
                if edl:
                    doc += "\n    Supported values:"
                    for ed in edl:
                        doc += "\n        " + str(ed.value) + " = '"+ed.name.decode("utf-8") + "'"
            return doc
        return "No documentation available."

    def __str__(self):
        _name = self.name
        if _name:
            return _name.decode("utf-8")
        else:
            return "<no name provided>"
        return None

    def __eq__(self, value):
        return self.__str__() == value.__str__()

    def __getattr__(self, name):
        field_descriptor = self.callbacks.field_descriptor_list.field_descriptor_for_name(name)
        return REB_BINARYDATA_DTYPE[field_descriptor.dtype].from_address(self.state + field_descriptor.offset).value

    def __setattr__(self, name, value):
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
        # TODO: implement proper output
        if self.name:
            return '<{0}.{1} object at {2}, name=\'{3}\'>'.format(self.__module__, type(self).__name__, hex(id(self)), self.name.decode("ascii"))
        else:
            return '<{0}.{1} object at {2}, name=None>'.format(self.__module__, type(self).__name__, hex(id(self)))

