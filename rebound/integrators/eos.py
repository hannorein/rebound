import ctypes


EOS_TYPES = {
        "lf": 0x00,
        "lf4": 0x01,
        "lf6": 0x02,
        "lf8": 0x03,
        "lf4_2": 0x04,
        "lf8_6_4": 0x05,
        "plf7_6_4": 0x06,
        "pmlf4": 0x07,
        "pmlf6": 0x08,
        }


class IntegratorEOS(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_eos.
    It controls the behaviour of the Embedded Operator Splitting methods. See Rein (2019) 
    for more details.
    
    :ivar int,string phi0      
        Sets the Phi_0 operator splitting method
    :ivar int,string phi1     
        Sets the Phi_1 operator splitting method
    :ivar int n     
        Sets the number of substeps taken by Phi_1
    :ivar int safe_mode  
        By default, safe_mode is on (1). Set to 0 (off) to combine
        drift step at the beginning and end of the Phi0 integrator steps.

    Example usage:
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "eos"
    >>> sim.ri_eos.phi0 = "LF8_6_4"
    >>> sim.ri_eos.phi1 = "LF8"
    >>> sim.ri_eos.n = 1
    >>> sim.ri_eos.safe_mode = 0

    """
    def __repr__(self):
        return '<{0}.{1} object at {2}, safe_mode={3}, is_synchronized={4}, n={5}, phi0={6}, phi2={7}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.safe_mode, self.is_synchronized, self.n, self.phi0, self.phi1)

    @property
    def phi0(self):
        """
        Get or set the type of operator splitting type for phi0.
        """
        i = self._phi0
        for name, _i in EOS_TYPES.items():
            if i==_i:
                return name
        return i
    @phi0.setter
    def phi0(self, value):
        if isinstance(value, int):
            self._phi0 = value
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "").replace("(", "").replace(")", "")
            if value in EOS_TYPES: 
                self._phi0 = EOS_TYPES[value]
            else:
                raise ValueError("Warning. EOS type %s not found."%value)
    @property
    def phi1(self):
        """
        Get or set the type of operator splitting type for phi1.
        """
        i = self._phi1
        for name, _i in EOS_TYPES.items():
            if i==_i:
                return name
        return i
    @phi1.setter
    def phi1(self, value):
        if isinstance(value, int):
            self._phi1 = value
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "").replace("(", "").replace(")", "")
            if value in EOS_TYPES: 
                self._phi1 = EOS_TYPES[value]
            else:
                raise ValueError("Warning. EOS type %s not found."%value)
    _fields_ = [
                ("_phi0",ctypes.c_uint),
                ("_phi1",ctypes.c_uint),
                ("n",ctypes.c_uint),
                ("safe_mode",ctypes.c_uint),
                ("is_synchronized",ctypes.c_uint),
                ]

