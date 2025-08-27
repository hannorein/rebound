from ctypes import POINTER, Structure

class Orbithierarchy(Structure):
    """
    A class representing one orbit in a tree. 
    """
    _fields_ = [("primary", POINTER(Orbithierarchy)),
                ("secondary", POINTER(Orbithierarchy)),
                ("com", POINTER(Particle)),
                ]

    def __repr__(self):
        return "<rebound.Orbithierarchy instance, primary={0} secondary={1} com={2}>".format(str(self.primary),str(self.secondary), str(self.com))


