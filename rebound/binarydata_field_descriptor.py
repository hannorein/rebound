import ctypes
from . import clibrebound, string_size_max
from .vectors import Vec3dBasic
import re
from textwrap import TextWrapper


# Note: ignoring some types that are not supposed to be set by the user
REB_BINARYDATA_DTYPE = {1: ctypes.c_double, 2: ctypes.c_int, 3: ctypes.c_uint, 4: ctypes.c_uint32, 5: ctypes.c_int64,
                        6: ctypes.c_uint64, 7: ctypes.c_size_t, 8: Vec3dBasic, 10: ctypes.c_void_p, 11: ctypes.c_void_p}

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

class DescriptorListIterator:
    """ Iterator for both BindarydataFieldDescriptor and EnumDescriptor lists """
    def __init__(self, dlist):
        self._dlist = dlist
        self._index = 0
    def __next__(self):
        if not ctypes.cast(self._dlist, ctypes.c_void_p).value:
            # Don't iterate NULL pointer
            raise StopIteration
        ed = self._dlist[self._index]
        if ed.name == b"":
            raise StopIteration
        self._index += 1 
        return ed
    def __iter__(self):
        return self

class EnumDescriptor(ctypes.Structure):
    """ Describes an Enum value in C """
    _fields_ = [("value", ctypes.c_int), 
                ("name", ctypes.c_char*string_size_max),
                ]
BaseEnumDescriptorList = ctypes.POINTER(EnumDescriptor)
class EnumDescriptorList(BaseEnumDescriptorList):
    _type_ = EnumDescriptor
    def __iter__(self):
        return DescriptorListIterator(self)


class BinarydataFieldDescriptor(ctypes.Structure):
    """ Describes the binarydata field in simulationarchives  """
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

    def doc(self):
        if self.name == b'':
            return None
        if self.documentation == b'':
            return None
        doc = self.name.decode("utf-8") + " : " 
        if self.dtype not in REB_BINARYDATA_DTYPE:
            doc += "(unknown datatype)\n"
        else:
            doc += str(REB_BINARYDATA_DTYPE[self.dtype].__name__) + "\n"
        doc += format_doc(self.documentation,indent="    ")
        edl = self.enum_descriptor_list
        if edl:
            doc += "\n    Supported values:"
            for ed in edl:
                doc += "\n        " + str(ed.value) + " = '"+ed.name.decode("utf-8") + "'"
        return doc

BaseBinarydataFieldDescriptorList = ctypes.POINTER(BinarydataFieldDescriptor)
class BinarydataFieldDescriptorList(BaseBinarydataFieldDescriptorList):
    _type_ = BinarydataFieldDescriptor
    def __iter__(self):
        return DescriptorListIterator(self)

    def field_descriptor_for_name(self, name):
        for field_descriptor in self: 
            if field_descriptor.name.decode("utf-8") == name:
                if field_descriptor.dtype not in REB_BINARYDATA_DTYPE:
                    raise NotImplemented("Datatype of field '%s' in IntegratorState is currently not supported." % name)
                return field_descriptor
        raise AttributeError("Field '%s' not found in BinarydataFieldDescriptorList." % name)
