import ctypes

SABA_TYPES = {
        "1": 0x0, "2": 0x1, "3": 0x2, "4": 0x3,
        "cm1": 0x100, "cm2": 0x101, "cm3": 0x102, "cm4": 0x103,
        "cl1": 0x200, "cl2": 0x201, "cl3": 0x202, "cl4": 0x203,
        "10,4": 0x4, "8,6,4": 0x5, "10,6,4": 0x6,
        "h8,4,4": 0x7, "h8,6,4": 0x8, "h10,6,4": 0x9,
        }

class IntegratorSABA(ctypes.Structure):
    """
    This class is an abstraction of the C-struct reb_integrator_saba.
    It controls the behaviour of the SABA integrator family.
    See Rein, Tamayo, and Brown (2019) for more details.

    :ivar str type:      
        Set the type of SABA integrator manually. The type can also be set by setting
        the integrator field in the REBOUND simulation. 

    :ivar int safe_mode:      
        This variable acts the same as for WHFast.
        If safe_mode is 1 (default) particles can be modified between
        timesteps and particle velocities and positions are always synchronised.
        If you set safe_mode to 0, the speed and accuracy of the integrator will improve.
        However, make sure you are aware of the consequences. Read the iPython tutorial
        on advanced WHFast usage to learn more.
   
    Example usage:
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "SABA(10,6,4)"
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "SABA"
    >>> sim.ri_saba.type = "(10,6,4)"
    
    >>> sim = rebound.Simulation()
    >>> sim.integrator = "SABACL4"
    >>> sim.ri_saba.safe_mode = 0

    """
    _fields_ = [("_type", ctypes.c_uint),
                ("safe_mode", ctypes.c_uint),
                ("is_synchronized", ctypes.c_uint),
                ("keep_unsynchronized", ctypes.c_uint),
            ]
    @property
    def type(self):
        """
        Get or set the type of SABA integrator.
        """
        i = self._type
        for name, _i in SABA_TYPES.items():
            if i==_i:
                return name
        return i
    @type.setter
    def type(self, value):
        if isinstance(value, int):
            self._type = value
        elif isinstance(value, basestring):
            value = value.lower().replace(" ", "").replace("(", "").replace(")", "")
            if value in SABA_TYPES: 
                self._type = SABA_TYPES[value]
            else:
                raise ValueError("Warning. SABA type not found.")

