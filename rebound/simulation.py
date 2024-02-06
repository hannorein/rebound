from ctypes import Structure, c_double, POINTER, c_uint32, c_float, c_int, c_uint, c_int64, c_uint64, c_void_p, c_char_p, CFUNCTYPE, byref, create_string_buffer, addressof, pointer, c_char, c_size_t, string_at 
from . import clibrebound, Escape, NoParticles, Encounter, Collision, GenericError 
from .citations import cite
from .units import units_convert_particle, check_units, convert_G, hash_to_unit
from .hash import hash as rebhash, HashPointerPair
from .vectors import Vec3d, Vec3dBasic, Vec6d
import os
import sys
import warnings

import types
      
### The following enum and class definitions need to
### consitent with those in rebound.h
        
INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "leapfrog": 4, "none": 7, "janus": 8, "mercurius": 9, "saba": 10, "eos": 11, "bs": 12, "whfast512":21}
BOUNDARIES = {"none": 0, "open": 1, "periodic": 2, "shear": 3}
GRAVITIES = {"none": 0, "basic": 1, "compensated": 2, "tree": 3, "mercurius": 4, "jacobi": 5}
COLLISIONS = {"none": 0, "direct": 1, "tree": 2, "line": 4, "linetree": 5}
# Format: Majorerror, id, message
BINARY_WARNINGS = [
    (True,  1, "Cannot read binary file. Check filename and file contents."),
    (False, 2, "Binary file was saved with a different version of REBOUND. Binary format might have changed."),
    (False, 4, "You have to reset function pointers after creating a reb_simulation struct with a binary file."),
    (False, 8, "Binary file might be corrupted. Number of particles found does not match particle number expected."),
    (True,  16, "Error while reading binary file (file was closed).",),
    (True,  32, "Index out of range.",),
    (True,  64, "Error while trying to seek file.",),
    (False, 128, "Encountered unkown field in file. File might have been saved with a different version of REBOUND."),
    (True,  256, "Integrator type is not supported by this simulationarchive version."),
    (False,  512, "The binary file seems to be corrupted. An attempt has been made to read the uncorrupted parts of it."),
    (True, 1024, "Reading old Simulationarchives (version < 2) is no longer supported. If you need to read such an archive, use a REBOUND version <= 3.26.3"),
]

# Note: name conflict with exception "Collision"
class CollisionS(Structure):
    _fields_ = [("p1", c_int),
                ("p2", c_int),
                ("gb", Vec6d),
                ("ri", c_int)]
    
    def __repr__(self):
        return '<{0}.{1} object at {2}, p1={3}, p2={4}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.p1, self.p2)
    


class Simulation(Structure):
    """
    This is the REBOUND Simulation Class.
    In encapsulated an entire REBOUND simulation and is an abstraction of the C struct reb_simulation.

    ### Examples
    
    Most simulation parameters can be directly changed with the property syntax:

    >>> sim = rebound.Simulation()
    >>> sim.G = 1.                  # Sets the graviational constant (default 1)
    >>> sim.softening = 1.          # Sets the graviational softening parameter (default 0)
    >>> sim.testparticle_type = 1   # Allows massive particles to feel influence from testparticles (default 0)
    >>> sim.dt = 0.1                # Sets the timestep (will change for adaptive integrators such as IAS15).
    >>> sim.t = 0.                  # Sets the current simulation time (default 0)
    >>> print(sim.N)                # Gets the current number of particles
    >>> print(sim.N_active)         # Gets the current number of active particles
       
    By calling rebound.Simulation() as shown above, you create a new simulation object
    The following example creates a simulation, saves it to a file and then creates
    a copy of the simulation stored in the binary file.

    >>> sim = rebound.Simulation()
    >>> sim.add(m=1.)
    >>> sim.add(m=1.e-3,x=1.,vy=1.)
    >>> sim.save_to_file("simulation.bin")
    >>> sim_copy = rebound.Simulation("simulation.bin")
    
    Similarly, you can create a simulation from a Simulationarchive
    by specifying the snapshot you want to load. 

    >>> sim = rebound.Simulation("archive.bin", 34)
    
    or 
    
    >>> sim = rebound.Simulation(filename="archive.bin", snapshot=34)

    Finally, you can also create a new Simulation by passing a bytes object of a 
    Simulationarchive to Simulation():
    
    >>> sim = rebound.Simulation(open("archive.bin","rb").read())

    """
    def __new__(cls, *args, **kw):
        # Create a new simulation if no arguments given
        if len(args)==0:
            sim = super(Simulation,cls).__new__(cls)
            clibrebound.reb_simulation_init(byref(sim))
            return sim
        
        # If first argument is of type bytes, then unpickle a Simulation
        if isinstance(args[0], bytes):
            l = len(args[0]) 
            buft = c_char * l
            buf = buft.from_buffer_copy(args[0])
            # Note: Not calling Simulationarchive.
            # Doing this manually here because we need to keep the reference to buf.
            # So we can later access the contents of the archive to get the simulation.
            sa = Simulationarchive.__new__(Simulationarchive, None, None)
            w = c_int(0)
            clibrebound.reb_simulationarchive_init_from_buffer_with_messages(byref(sa), byref(buf), c_size_t(l), None, byref(w))
            sim = super(Simulation,cls).__new__(cls)
            clibrebound.reb_simulation_init(byref(sim))
            clibrebound.reb_simulation_create_from_simulationarchive_with_messages(byref(sim),byref(sa),c_int64(-1),byref(w))
            for majorerror, value, message in BINARY_WARNINGS:
                if w.value & value:
                    if majorerror:
                        raise RuntimeError(message)
                    else:  
                        # Just a warning
                        warnings.warn(message, RuntimeWarning)
            return sim
       
        # Create simulation from Simulationarchive
        if isinstance(args[0], Simulationarchive):
            sa = args[0]
        else:
            # Otherwise assume first argument is filename
            filename = args[0]
            if "filename" in kw:
                filename = kw["filename"]
            sa = Simulationarchive(filename,process_warnings=False)

        snapshot = -1
        if len(args)>1:
            snapshot = args[1]
        if "snapshot" in kw:
            snapshot = kw["snapshot"]

        if sa is not None:
            # Recreate exisitng simulation 
            sim = super(Simulation,cls).__new__(cls)
            clibrebound.reb_simulation_init(byref(sim))
            w = sa.warnings # warnings will be appended to previous warnings (as to not repeat them) 
            clibrebound.reb_simulation_create_from_simulationarchive_with_messages(byref(sim),byref(sa),c_int64(snapshot),byref(w))
            for majorerror, value, message in BINARY_WARNINGS:
                if w.value & value:
                    if majorerror:
                        raise RuntimeError(message)
                    else:  
                        # Just a warning
                        warnings.warn(message, RuntimeWarning)
            return sim

        # Still here? Then an error occured.
        raise RuntimeError("Can not create Simulation.")

    def __init__(self,filename=None,snapshot=None):
        self.save_messages = 1 # Warnings will be checked within python

    def __repr__(self):
        return '<{0}.{1} object at {2}, N={3}, t={4}>'.format(self.__module__, type(self).__name__, hex(id(self)), self.N, self.t)
    
    @classmethod
    def from_simulationarchive(cls, simulationarchive, snapshot=-1):
        return cls(filename=filename,snapshot=snapshot)

    @classmethod
    def from_file(cls, filename, snapshot=-1):
        return cls(filename=filename, snapshot=snapshot)
    
    def copy(self):
        """
        Returns a deep copy of a REBOUND simulation. You need to reset 
        any function pointers on the copy. 
        
        Returns
        ------- 
        A rebound.Simulation object.
        
        """
        w = c_int(0)
        sim = Simulation()
        clibrebound.reb_simulation_copy_with_messages(byref(sim),byref(self),byref(w))
        for majorerror, value, message in BINARY_WARNINGS:
            if w.value & value:
                if majorerror:
                    raise RuntimeError(message)
                else:  
                    # Just a warning
                    warnings.warn(message, RuntimeWarning)
        return sim

    def start_server(self, port=1234):
        """
        Start webserver on specified port. 
        The default port is 1234.
        You can access a running server by opening a web browser
        at http://localhost:1234 or http://127.0.0.1:1234
        """
        clibrebound.reb_simulation_start_server.restype = c_int
        ret_value = clibrebound.reb_simulation_start_server(byref(self), c_int(port))
        self.process_messages()
    
    def stop_server(self, port=1234):
        """
        Stop the webserver.
        """
        ret_value = clibrebound.reb_simulation_stop_server(byref(self))
        self.process_messages()


    def widget(self, port=None, host="localhost", size=(500,500)):
        if port is not None:
            port = int(port)
            if self._server_data:
                if self._server_data.contents.port != port:
                    raise RuntimeError("Server is running already on port %d but port %d was requested."%(self._server_data.contents.port, port))
        if not self._server_data:
            if port is not None:
                self.start_server(port=port)
            else:
                self.start_server()
        if not self._server_data:
            raise RuntimeError("Server did not start up immediately. Try again in a few seconds.")
        port = int(self._server_data.contents.port)
        from IPython.display import IFrame
        width, height = size
        display(IFrame("http://"+host+":"+"%d"%port, width, height))

    def cite(self):
        """
        Generate citations

        This function generates citations to papers relevant to the current 
        setting of the simulation.
        """

        txt, bib = cite(self)
        # one could check for REBOUNDx here, then append txt and bib accordingly
        print(txt + "\n\n\n" + bib)

    @property
    def simulationarchive_filename(self):
        """
        Returns the current Simulationarchive filename in use. 
        Do not set manually. Use sim.save_to_file() instead
        """
        return self._simulationarchive_filename

# Message and memory management functions
    def process_messages(self):
        clibrebound.reb_simulation_get_next_message.restype = c_int
        buf = create_string_buffer(c_int.in_dll(clibrebound, "reb_max_messages_length").value)
        while clibrebound.reb_simulation_get_next_message(byref(self), buf):
            msg = buf.value.decode("ascii")
            if msg[0]=='w':
                warnings.warn(msg[1:], RuntimeWarning)
            elif msg[0]=='e':
                raise RuntimeError(msg[1:])

# Pickling methods: return Simulationarchive binary
    def __reduce__(self):
        buf = c_char_p()
        size = c_size_t()
        clibrebound.reb_simulation_save_to_stream(byref(self), byref(buf), byref(size))
        s = bytes(string_at(buf, size=size.value)) #make copy
        clibrebound.reb_simulation_output_free_stream(buf) # free original
        return (Simulation, (s,))

# Other operators

    def __del__(self):
        if self._b_needsfree_ == 1: # to avoid, e.g., sim.particles[1]._sim.contents.G creating a Simulation instance to get G, and then freeing the C simulation when it immediately goes out of scope
            clibrebound.reb_simulation_free_pointers(byref(self))

    def __eq__(self, other):
        # This ignores the walltime parameter
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_simulation_diff.restype = c_int
        ret = clibrebound.reb_simulation_diff(byref(self), byref(other),c_int(2))
        return not ret
            
    def diff(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_simulation_diff_char.restype = c_char_p
        output = clibrebound.reb_simulation_diff_char(byref(other), byref(self))
        print(output.decode("utf-8"))

    def __add__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        c = self.copy()
        return c.__iadd__(other)
    
    def __iadd__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_simulation_iadd.restype = c_int
        ret = clibrebound.reb_simulation_iadd(byref(self), byref(other))
        if ret==-1:
            raise RuntimeError("Cannot add simulations. Check that the simulations have the same number of particles")
        return self
    
    def __sub__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        c = self.copy()
        return c.__isub__(other)

    def __isub__(self, other):
        if not isinstance(other,Simulation):
            return NotImplemented
        clibrebound.reb_simulation_isub.restype = c_int
        ret = clibrebound.reb_simulation_isub(byref(self), byref(other))
        if ret==-1:
            raise RuntimeError("Cannot subtract simulations. Check that the simulations have the same number of particles")
        return self
    
    def __mul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        c = self.copy()
        c.multiply(other, other)
        return c
    
    def __imul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        self.multiply(other, other)
        return self

    def __rmul__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        c = self.copy()
        c.multiply(other, other)
        return c
    
    def __div__(self, other):
        return self.__truediv__(other)
    
    def __idiv__(self, other):
        return self.__itruediv__(other)

    def __truediv__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        c = self.copy()
        if other==0.:
            raise ZeroDivisionError
        c.multiply(1./other, 1./other)
        return c

    def __itruediv__(self, other):
        try:
            other = float(other)
        except:
            return NotImplemented
        if other==0.:
            raise ZeroDivisionError
        self.multiply(1./other, 1./other)
        return self

    def multiply(self, scalar_pos, scalar_vel):
        try:
            scalar_pos = float(scalar_pos)
            scalar_vel = float(scalar_vel)
        except:
            raise ValueError("Cannot multiply simulation with non-scalars.")
        clibrebound.reb_simulation_imul(byref(self), c_double(scalar_pos), c_double(scalar_vel))
    
    def rotate(self, q):
        from .rotation import Rotation
        if not isinstance(q, Rotation):
            return NotImplemented
        clibrebound.reb_simulation_irotate(byref(self), q)

#ODE functions
    def create_ode(self, length, needs_nbody=True):
        clibrebound.reb_ode_create.restype = POINTER(ODE)
        ode_p = clibrebound.reb_ode_create(byref(self), c_int(length))
        ode_p.contents.needs_nbody = c_uint(needs_nbody)
        return ODE.from_address(addressof(ode_p.contents))

# Status functions
    def status(self, showParticles=True, showAllFields=True):
        """ 
        Prints a summary of the current status 
        of the simulation.
        """
        from rebound import __version__, __build__
        s= ""
        s += "---------------------------------\n"
        s += "REBOUND version:     \t%s\n" %__version__
        s += "REBOUND built on:    \t%s\n" %__build__
        s += "Number of particles: \t%d\n" %self.N       
        s += "Selected integrator: \t" + self.integrator + "\n"       
        s += "Simulation time:     \t%.16e\n" %self.t
        s += "Current timestep:    \t%f\n" %self.dt
        s += "---------------------------------\n"
        if self.N>0 and showParticles:
            for p in self.particles:
                s += str(p) + "\n"
            s += "---------------------------------\n"
        print(s, end="")
        if showAllFields:
            print("The following fields have non-default values:")
            newsim = Simulation()
            clibrebound.reb_simulation_diff_char.restype = c_char_p
            output = clibrebound.reb_simulation_diff_char(byref(newsim), byref(self))
            print(output.decode("utf-8"))



# Set function pointer for additional forces
    @property
    def additional_forces(self):
        """
        Get or set a function pointer for calculating additional forces in the simulation.

        The argument can be a python function or something that can 
        be cast to a C function of type CFUNCTYPE(None,POINTER(Simulaton)). 
        If the forces are velocity dependent, the flag 
        force_is_velocity_dependent needs to be set to 1. Otherwise,
        the particle structures might contain incorrect velocity 
        values.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @additional_forces.setter
    def additional_forces(self, func):
        self._afp = AFF(func)
        self._additional_forces = self._afp

    @property
    def pre_timestep_modifications(self):
        """
        Get or set a function pointer for pre-timestep modifications.

        The argument can be a python function or something that can be cast to a C function or a
        python function.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @pre_timestep_modifications.setter
    def pre_timestep_modifications(self, func):
        self._pretmp = AFF(func)
        self._pre_timestep_modifications = self._pretmp
    
    @property
    def post_timestep_modifications(self):
        """
        Get or set a function pointer for post-timestep modifications.

        The argument can be a python function or something that can be cast to a C function or a
        python function.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @post_timestep_modifications.setter
    def post_timestep_modifications(self, func):
        self._posttmp = AFF(func)
        self._post_timestep_modifications = self._posttmp
 
    @property
    def heartbeat(self):
        """
        Set a function pointer for a heartbeat function.
        The heartbeat function is called every timestep and can be used
        to monitor long simulations, check for stalled simulations and 
        output debugging information.
     
        The argument can be a python function or something that can be 
        cast to a C function or a python function.

        The function called will receive a pointer to the simulation
        object as its argument.
        
        Examples
        --------

        >>> import rebound
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(m=1e-3, a=1)
        >>> def heartbeat(sim):
        >>>     # sim is a pointer to the simulation object,
        >>>     # thus use contents to access object data.
        >>>     # See ctypes documentation for details.
        >>>     print(sim.contents.t)
        >>> sim.heartbeat = heartbeat
        >>> sim.integrate(1.)

        """
        raise AttributeError("You can only set C function pointers from python (not get).")
    @heartbeat.setter
    def heartbeat(self, func):
        self._hb = AFF(func)
        self._heartbeat = self._hb

    @property 
    def coefficient_of_restitution(self):
        """
        Get or set a function pointer that defined the coefficient of restitution.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @coefficient_of_restitution.setter
    def coefficient_of_restitution(self, func):
        self._corfp = CORFF(func)
        self._coefficient_of_restitution = self._corfp
    
    @property 
    def collision_resolve(self):
        """
        Get or set a function pointer for collision resolving routine.
        
        Possible options for setting:
          1) Function pointer
          2) "merge": two colliding particles will merge) 
          3) "hardsphere": two colliding particles will bounce off using a set coefficient of restitution
        """
        raise AttributeError("You can only set C function pointers from python.")
    @collision_resolve.setter
    def collision_resolve(self, func):
        if func == "merge":
            clibrebound.reb_simulation_set_collision_resolve.restype = None
            clibrebound.reb_simulation_set_collision_resolve(byref(self), clibrebound.reb_collision_resolve_merge)
        elif func == "hardsphere":
            clibrebound.reb_simulation_set_collision_resolve.restype = None
            clibrebound.reb_simulation_set_collision_resolve(byref(self), clibrebound.reb_collision_resolve_hardsphere)
        elif func == "halt":
            clibrebound.reb_simulation_set_collision_resolve.restype = None
            clibrebound.reb_simulation_set_collision_resolve(byref(self), clibrebound.reb_collision_resolve_halt)
        else:
            self._colrfp = COLRFF(func)
            self._collision_resolve = self._colrfp
    
    @property 
    def free_particle_ap(self):
        """
        Get or set a function pointer for freeing the ap pointer whenever sim.remove is called.
        """
        raise AttributeError("You can only set C function pointers from python.")
    @free_particle_ap.setter
    def free_particle_ap(self, func):
        self._fpa = FPA(func)
        self._free_particle_ap = self._fpa

# Setter/getter of parameters and constants
    @property 
    def N_real(self):
        """
        Get the current number of real (i.e. no variational/shadow) particles in the simulation.
        """
        return self.N-self.N_var

    @property
    def integrator(self):
        """
        Get or set the integrator module.

        Available integrators include:

        - ``'IAS15'`` (default)
        - ``'WHFast'``
        - ``'SEI'``
        - ``'LEAPFROG'``
        - ``'JANUS'``
        - ``'MERCURIUS'``
        - ``'WHCKL'`` 
        - ``'WHCKM'`` 
        - ``'WHCKC'`` 
        - ``'SABA4'`` 
        - ``'SABACL4'`` 
        - ``'SABACM4'`` 
        - ``'SABA(10,6,4)'`` 
        - ``'EOS'`` 
        - ``'BS'`` 
        - ``'WHFast512'``
        - ``'none'``
        
        Check the online documentation for a full description of each of the integrators. 
        """
        i = self._integrator
        for name, _i in INTEGRATORS.items():
            if i==_i:
                return name
        return i
    @integrator.setter
    def integrator(self, value):
        if isinstance(value, int):
            self._integrator = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in INTEGRATORS: 
                self._integrator = INTEGRATORS[value]
            # Shortcuts
            elif value=="wh":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 0
                self.ri_whfast.kernel = "default"
            elif value=="whc":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "default"
            elif value=="whckl":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "lazy"
            elif value=="whckm":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "modifiedkick"
            elif value=="whckc":
                self.integrator = "whfast"
                self.ri_whfast.corrector = 17
                self.ri_whfast.kernel = "composition"
            elif value[0:4]=="saba" and len(value)>4:
                self.integrator = "saba"
                self.ri_saba.type = value[4:]
            else:
                raise ValueError("Integrator not found.")
    
    @property
    def boundary(self):
        """
        Get or set the boundary module.

        Available boundary modules are:

        - ``'none'`` (default)
        - ``'open'``
        - ``'periodic'``
        - ``'shear'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._boundary
        for name, _i in BOUNDARIES.items():
            if i==_i:
                return name
        return i
    @boundary.setter
    def boundary(self, value):
        if isinstance(value, int):
            self._boundary = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in BOUNDARIES: 
                self._boundary = BOUNDARIES[value]
            else:
                raise ValueError("Warning. Boundary condition module not found.")

    @property
    def gravity(self):
        """
        Get or set the gravity module.

        Available gravity modules are:

        - ``'none'`` 
        - ``'basic'`` (default)
        - ``'compensated'``
        - ``'tree'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._gravity
        for name, _i in GRAVITIES.items():
            if i==_i:
                return name
        return i
    @gravity.setter
    def gravity(self, value):
        if isinstance(value, int):
            self._gravity = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in GRAVITIES: 
                self._gravity = GRAVITIES[value]
            else:
                raise ValueError("Warning. Gravity module not found.")

    @property
    def collision(self):
        """
        Get or set the collision module.

        Available collision modules are:

        - ``'none'`` (default)
        - ``'direct'``
        - ``'tree'``
        - ``'line'`` 
        - ``'linetree'``
        
        Check the online documentation for a full description of each of the modules. 
        """
        i = self._collision
        for name, _i in COLLISIONS.items():
            if i==_i:
                return name
        return i
    @collision.setter
    def collision(self, value):
        if isinstance(value, int):
            self._collision = c_int(value)
        elif isinstance(value, basestring):
            value = value.lower()
            if value in COLLISIONS: 
                self.collision = COLLISIONS[value]
            else:
                raise ValueError("Warning. Collision module not found.")

# Units

    @property
    def units(self):
        """
        Tuple of the units for length, time and mass.  Can be set in any order, and strings are not case-sensitive.  See ipython_examples/Units.ipynb for more information.  You can check the units' exact values and add Additional units in rebound/rebound/units.py.  Units should be set before adding particles to the simulation (will give error otherwise).

        Currently supported Units 
        -------------------------

        Times:
        Hr          : Hours
        Yr          : Julian years
        Jyr         : Julian years
        Sidereal_yr : Sidereal year
        Yr2pi       : Year divided by 2pi, with year defined as orbital period of planet at 1AU around 1Msun star
        Kyr         : Kiloyears (Julian)
        Myr         : Megayears (Julian)
        Gyr         : Gigayears (Julian)

        Lengths:
        M           : Meters
        Cm          : Centimeters
        Km          : Kilometers
        AU          : Astronomical Units

        Masses:
        Kg          : Kilograms
        Msun        : Solar masses
        Mmercury    : Mercury masses
        Mvenus      : Venus masses
        Mearth      : Earth masses
        Mmars       : Mars masses
        Mjupiter    : Jupiter masses
        Msaturn     : Saturn masses
        Muranus     : Neptune masses
        Mpluto      : Pluto masses
        
        Examples
        --------
        
        >>> sim = rebound.Simulation()
        >>> sim.units = ('yr', 'AU', 'Msun')

        """
        
        return {'length':hash_to_unit(self.python_unit_l), 'mass':hash_to_unit(self.python_unit_m), 'time':hash_to_unit(self.python_unit_t)}

    @units.setter
    def units(self, newunits):
        newunits = check_units(newunits)        
        if self.N>0: # some particles are loaded
            raise AttributeError("Error:  You cannot set the units after populating the particles array.  See ipython_examples/Units.ipynb.")
        self.update_units(newunits) 

    def update_units(self, newunits): 
        clibrebound.reb_hash.restype = c_uint32
        self.python_unit_l = clibrebound.reb_hash(c_char_p(newunits[0].encode("ascii")))
        self.python_unit_t = clibrebound.reb_hash(c_char_p(newunits[1].encode("ascii")))
        self.python_unit_m = clibrebound.reb_hash(c_char_p(newunits[2].encode("ascii")))
        self.G = convert_G(newunits)

    def convert_particle_units(self, *args): 
        """
        Will convert the units for the simulation (i.e. convert G) as well as the particles' cartesian elements.
        Must have set sim.units ahead of calling this function, so REBOUND knows what units to convert from.

        Parameters
        ----------
        3 strings corresponding to units of time, length and mass.  Can be in any order and aren't case-sensitive.        You can add new units to rebound/rebound/units.py
        """
        if self.python_unit_l == 0 or self.python_unit_m == 0 or self.python_unit_t == 0:
            raise AttributeError("Must set sim.units before calling convert_particle_units in order to know what units to convert from.")
        new_l, new_t, new_m = check_units(args)
        for p in self.particles:
            units_convert_particle(p, hash_to_unit(self.python_unit_l), hash_to_unit(self.python_unit_t), hash_to_unit(self.python_unit_m), new_l, new_t, new_m)
        self.update_units((new_l, new_t, new_m))

# Variational Equations
    def add_variation(self,order=1,first_order=None, first_order_2=None, testparticle=-1):
        """ 
        This function adds a set of variational particles to the simulation. 

        If there are N real particles in the simulation, this functions adds N additional variational 
        particles. To see how many particles (real and variational) are in a simulation, use ``'sim.N'``. 
        To see how many variational particles are in a simulation use ``'sim.N_var'``.

        Currently Leapfrog, WHFast and IAS15 support first order variational equations. IAS15 also
        supports second order variational equations.

        Parameters
        ----------
        order : integer, optional
            By default, the function adds a set of first order variational particles to the simulation. Set this flag to 2 for second order.
        first_order : Variation, optional
            Second order variational equations depend on their corresponding first order variational equations. 
            This parameter expects the Variation object corresponding to the first order variational equations. 
        first_order_2 : Variation, optional
            Same as first_order. But allows to set two different indicies to calculated off-diagonal elements. 
            If omitted, then first_order will be used for both first order equations.
        testparticle : int, optional
            If set to a value >= 0, then only one variational particle will be added and be treated as a test particle.
            

        Returns
        -------
        Returns Variation object (a copy--you can only modify it through its particles property or vary method). 
        """
        cur_N_var_config = self.N_var_config
        if order==1:
            index = clibrebound.reb_simulation_add_variation_1st_order(byref(self),c_int(testparticle))
        elif order==2:
            if first_order is None:
                raise AttributeError("Please specify corresponding first order variational equations when initializing second order variational equations.")
            if first_order_2 is None:
                first_order_2 = first_order
            index = clibrebound.reb_simulation_add_variation_2nd_order(byref(self),c_int(testparticle),c_int(first_order.index),c_int(first_order_2.index))
        else:
            raise AttributeError("Only variational equations of first and second order are supported.")

        # Need a copy because location of original might shift if more variations added
        s = Variation.from_buffer_copy(self.var_config[cur_N_var_config])

        return s
        
# MEGNO
    def init_megno(self, seed=None):
        """
        This function initialises the chaos indicator MEGNO particles and enables their integration.

        MEGNO is short for Mean Exponential Growth of Nearby orbits. It can be used to test
        if a system is chaotic or not. In the backend, the integrator is integrating an additional set
        of particles using the variational equation. Note that variational equations are better 
        suited for this than shadow particles. MEGNO is currently only supported in the IAS15 
        and WHFast integrators.

        This function also needs to be called if you are interested in the Lyapunov exponent as it is
        calculate with the help of MEGNO. See Rein and Tamayo 2015 for details on the implementation.

        For more information on MEGNO see e.g. https://dx.doi.org/10.1051/0004-6361:20011189
        """
        if seed is None:
            clibrebound.reb_simulation_init_megno(byref(self))
        else:
            clibrebound.reb_simulation_init_megno_seed(byref(self), c_uint(seed))
    
    def megno(self):
        """
        Return the current MEGNO value.
        Note that you need to call init_megno() before the start of the simulation.
        """
        if self._calculate_megno==0:
            raise RuntimeError("MEGNO cannot be calculated. Make sure to call init_megno() after adding all particles but before integrating the simulation.")

        clibrebound.reb_simulation_megno.restype = c_double
        return clibrebound.reb_simulation_megno(byref(self))
    
    def lyapunov(self):
        """
        Return the current Lyapunov Characteristic Number (LCN).
        Note that you need to call init_megno() before the start of the simulation.
        Different definitions of the LCN exist.  Here, we're following Eq 24 of 
        Cincotta and Simo (2000): https://aas.aanda.org/articles/aas/abs/2000/20/h1686/h1686.html.
        To get a timescale (the Lyapunov timescale), take the inverse of this quantity.
        """
        if self._calculate_megno==0:
            raise RuntimeError("Lyapunov Characteristic Number cannot be calculated. Make sure to call init_megno() after adding all particles but before integrating the simulation.")

        clibrebound.reb_simulation_lyapunov.restype = c_double
        return clibrebound.reb_simulation_lyapunov(byref(self))
    
# Particle add function, used to be called particle_add() and add_particle() 
    def add(self, particle=None, **kwargs):   
        """
        Adds a particle to REBOUND. Accepts one of the following:

        1) A single Particle structure.
        2) The particle's mass and a set of cartesian coordinates: m,x,y,z,vx,vy,vz.
        3) The primary as a Particle structure, the particle's mass and a set of orbital elements: primary,m,a,anom,e,omega,inv,Omega,MEAN (see :class:`.Orbit` for the definition of orbital elements).
        4) A name of an object (uses NASA Horizons to look up coordinates)
        5) A list of particles or names.
        """
        if particle is not None:
            if isinstance(particle, Particle):
                if (self.gravity == "tree" or self.collision == "tree") and self.root_size <=0.:
                    raise ValueError("The tree code for gravity and/or collision detection has been selected. However, the simulation box has not been configured yet. You cannot add particles until the the simulation box has a finite size.")

                clibrebound.reb_simulation_add(byref(self), particle)
            elif isinstance(particle, list):
                for p in particle:
                    self.add(p, **kwargs)
            elif isinstance(particle,str):
                if self.python_unit_l == 0 or self.python_unit_m == 0 or self.python_unit_t == 0:
                    self.units = ('AU', 'yr2pi', 'Msun')
                    self.G = 1.0
                builtindatasets = ["solar system", "outer solar system"]
                if particle.lower() == "solar system":          # built in test dataset
                    data.add_solar_system(self)
                elif particle.lower() == "outer solar system":  # built in test dataset
                    data.add_outer_solar_system(self)
                else:
                    if "frame" not in kwargs:
                        if hasattr(self, 'default_plane'):
                            kwargs["plane"] = self.default_plane # allow ASSIST to set default plane
                    self.add(horizons.query_horizons_for_particle(particle, **kwargs), hash=particle)
                    units_convert_particle(self.particles[-1], 'km', 's', 'kg', hash_to_unit(self.python_unit_l), hash_to_unit(self.python_unit_t), hash_to_unit(self.python_unit_m))
            else: 
                raise ValueError("Argument passed to add() not supported.")
        else: 
            self.add(Particle(simulation=self, **kwargs))

# Particle getter functions
    @property
    def particles(self):
        """
        Returns a Particles object that allows users to access particles like a dictionary using indices, hashes, or strings. 

        The Particles object uses pointers and thus the contents update 
        as the simulation progresses. Note that the pointers could change,
        for example when a particle is added or removed from the simulation. 
        """
        particles = Particles(self)
        return particles

    @particles.deleter
    def particles(self):
        """
        Remove all particles from the simulation
        """
        clibrebound.reb_simulation_remove_all_particles(byref(self))
        self.process_messages()

    def remove(self, index=None, hash=None, keep_sorted=True):
        """ 
        Removes a particle from the simulation.

        Parameters
        ----------
        index : int, optional
            Specify particle to remove by index.
        hash : c_uint32 or string, optional
            Specifiy particle to remove by hash (if a string is passed, the corresponding hash is calculated).
        keep_sorted : bool, optional
            By default, remove preserves the order of particles in the particles array. 
            Might set it to zero in cases with many particles and many removals to speed things up.
        """
        if index is not None:
            clibrebound.reb_simulation_remove_particle(byref(self), index, keep_sorted)
        if hash is not None:
            hash_types = c_uint32, c_uint, c_uint64
            PY3 = sys.version_info[0] == 3
            if PY3:
                string_types = str,
                int_types = int,
            else:
                string_types = basestring,
                int_types = int, long
            if isinstance(hash, string_types):
                clibrebound.reb_simulation_remove_particle_by_hash(byref(self), rebhash(hash), keep_sorted)
            elif isinstance(hash, int_types):
                clibrebound.reb_simulation_remove_particle_by_hash(byref(self), c_uint32(hash), keep_sorted)
            elif isinstance(hash, hash_types):
                clibrebound.reb_simulation_remove_particle_by_hash(byref(self), hash, keep_sorted)

        self.process_messages()

# Orbit calculation
    def orbits(self, primary=None, jacobi_masses=False):
        """ 
        Calculate orbital parameters for all particles in the simulation.
        By default, this functions returns the orbits in Jacobi coordinates. 

        If MEGNO is enabled, variational particles will be ignored.

        Parameters
        ----------

        primary     : rebound.Particle, optional
            Set the primary against which to reference the osculating orbit. Default (use Jacobi center of mass).
            For heliocentric coordinates, pass the central object to this parameter. 
        jacobi_masses: bool
            Whether to use jacobi primary mass in orbit calculation. (Default: False)

        Returns
        -------
        Returns an array of Orbits of length N-1.
        """
        orbits = []
       
        if primary is None:
            jacobi = True
            primary = self.particles[0]
            clibrebound.reb_particle_com_of_pair.restype = Particle
        else:
            jacobi = False

        for p in self.particles[1:self.N_real]:
            if jacobi_masses is True:
                interior_mass = primary.m
                # orbit conversion uses mu=G*(p.m+primary.m) so set prim.m=Mjac-m so mu=G*Mjac
                primary.m = self.particles[0].m*(p.m + interior_mass)/interior_mass - p.m
                orbits.append(p.orbit(primary=primary))
                primary.m = interior_mass # back to total mass of interior bodies to update com
            else:
                orbits.append(p.orbit(primary=primary))
            if jacobi is True: # update com to include current particle for next iteration
                primary = clibrebound.reb_particle_com_of_pair(primary, p)

        return orbits

# COM calculation 
    def com(self, first=0, last=None):
        """
        Returns the center of momentum for all particles in the simulation.

        Parameters
        ----------
        first: int, optional
            If ``first`` is specified, only calculate the center of momentum starting
            from index=``first``.
        last : int or None, optional
            If ``last`` is specified only calculate the center of momentum up to 
            (but excluding) index=``last``.  Same behavior as Python's range function.

        Examples
        --------
        >>> sim = rebound.Simulation()
        >>> sim.add(m=1, x=-20)
        >>> sim.add(m=1, x=-10)
        >>> sim.add(m=1, x=0)
        >>> sim.add(m=1, x=10)
        >>> sim.add(m=1, x=20)
        >>> com = sim.com()
        >>> com.x
        0.0 
        >>> com = sim.com(first=2,last=4) # Considers indices 2,3
        >>> com.x
        5.0

        """
        if last is None:
            last = self.N_real
        clibrebound.reb_simulation_com_range.restype = Particle
        return clibrebound.reb_simulation_com_range(byref(self), c_int(first), c_int(last))

# Tools
    def serialize_particle_data(self,**kwargs):
        """
        Fast way to access serialized particle data via numpy arrays.

        This function can directly set the values of numpy arrays to
        current particle data. This is significantly faster than accessing
        particle data via `sim.particles` as all the copying is done 
        on the C side. 
        No memory is allocated by this function.
        It expects correctly sized numpy arrays as arguments. The argument
        name indicates what kind of particle data is written to the array. 
        
        Possible argument names are "hash", "m", "r", "xyz", "vxvyvz", and 
        "xyzvxvyvz". The datatype for the "hash" array needs to be uint32. 
        The other arrays expect a datatype of float64. The lengths of 
        "hash", "m", "r" arrays need to be at least sim.N. The lengths of 
        xyz and vxvyvz need to be at least 3*sim.N. The length of
        "xyzvxvyvz" arrays need to be 6*sim.N. Exceptions are raised 
        otherwise.

        Note that this routine is only intended for special use cases
        where speed is an issue. For normal use, it is recommended to
        access particle data via the `sim.particles` array. Be aware of
        potential issues that arise by directly accesing the memory of
        numpy arrays (see numpy documentation for more details).

        Examples
        --------
        This sets an array to the xyz positions of all particles:

        >>> import numpy as np
        >>> a = np.zeros((sim.N,3),dtype="float64")
        >>> sim.serialize_particle_data(xyz=a)
        >>> print(a)

        To get all current radii of particles:

        >>> a = np.zeros(sim.N,dtype="float64")
        >>> sim.serialize_particle_data(r=a)
        >>> print(a)
        
        To get all current radii and hashes of particles:

        >>> a = np.zeros(sim.N,dtype="float64")
        >>> b = np.zeros(sim.N,dtype="uint32")
        >>> sim.serialize_particle_data(r=a,hash=b)
        >>> print(a,b)

        """
        N = self.N
        possible_keys = ["hash","m","r","xyz","vxvyvz","xyzvxvyvz"]
        d = {x:None for x in possible_keys}
        for k,v in kwargs.items():
            if k in d:
                if k == "hash":
                    if v.dtype!= "uint32":
                        raise AttributeError("Expected 'uint32' data type for '%s' array."%k)
                    if v.size<N:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(POINTER(c_uint32))
                else:
                    if v.dtype!= "float64":
                        raise AttributeError("Expected 'float64' data type for %s array."%k)
                    if k in ["xyz", "vxvyvz"]:
                        minsize = 3*N
                    elif k in ["xyzvxvyvz"]:
                        minsize = 6*N
                    else:
                        minsize = N
                    if v.size<minsize:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(POINTER(c_double))
            else:
                raise AttributeError("Only '%s' are currently supported attributes for serialization." % "', '".join(d.keys()))

        clibrebound.reb_simulation_get_serialized_particle_data(byref(self), d["hash"], d["m"], d["r"], d["xyz"], d["vxvyvz"], d["xyzvxvyvz"])
    
    def set_serialized_particle_data(self,**kwargs):
        """
        Fast way to set serialized particle data via numpy arrays.
        This is the inverse of Simulation.serialize_particle_data()
        and uses the same syntax
        """

        N = self.N
        possible_keys = ["hash","m","r","xyz","vxvyvz","xyzvxvyvz"]
        d = {x:None for x in possible_keys}
        for k,v in kwargs.items():
            if k in d:
                if k == "hash":
                    if v.dtype!= "uint32":
                        raise AttributeError("Expected 'uint32' data type for '%s' array."%k)
                    if v.size<N:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(POINTER(c_uint32))
                else:
                    if v.dtype!= "float64":
                        raise AttributeError("Expected 'float64' data type for %s array."%k)
                    if k in ["xyz", "vxvyvz"]:
                        minsize = 3*N
                    elif k in ["xyzvxvyvz"]:
                        minsize = 6*N
                    else:
                        minsize = N
                    if v.size<minsize:
                        raise AttributeError("Array '%s' is not large enough."%k)
                    d[k] = v.ctypes.data_as(POINTER(c_double))
            else:
                raise AttributeError("Only '%s' are currently supported attributes for serialization." % "', '".join(d.keys()))

        clibrebound.reb_simulation_set_serialized_particle_data(byref(self), d["hash"], d["m"], d["r"], d["xyz"], d["vxvyvz"], d["xyzvxvyvz"])

    def move_to_hel(self):
        """
        This function moves all particles in the simulation to the heliocentric frame.
        Note that the simulation will not stay in the heliocentric frame during integrations
        as the heliocentric frame is not an inertial frame.
        """
        clibrebound.reb_simulation_move_to_hel(byref(self))
    
    def move_to_com(self):
        """
        This function moves all particles in the simulation to a center of momentum frame.
        In that frame, the center of mass is at the origin and does not move.
        It makes sense to call this function at the beginning of the integration, especially 
        for the high accuracy integrators IAS15 and WHFast.
        """
        clibrebound.reb_simulation_move_to_com(byref(self))

    def energy(self):
        """
        Returns the sum of potential and kinetic energy of all particles in the simulation.
        """
        clibrebound.reb_simulation_energy.restype = c_double
        return clibrebound.reb_simulation_energy(byref(self))
   
    def angular_momentum(self):
        """
        Returns a list of the three (x,y,z) components of the total angular momentum of all particles in the simulation.
        """
        clibrebound.reb_simulation_angular_momentum.restype = Vec3dBasic
        return Vec3d(clibrebound.reb_simulation_angular_momentum(byref(self)))

    def configure_box(self, boxsize, N_root_x=1, N_root_y=1, N_root_z=1):
        """
        Initialize the simulation box.

        This function only needs to be called it boundary conditions other than "none"
        are used. In such a case the boxsize must be known and is set with this function.

        Parameters
        ----------
        boxsize : float
            The size of one root box.
        N_root_x, N_root_y, N_root_z : int, optional
            The number of root boxes in each direction. The total size of the simulation box
            will be ``N_root_x * boxsize``, ``N_root_y * boxsize`` and ``N_root_z * boxsize``.
            By default, there will be exactly one root box in each direction.
        """
        clibrebound.reb_simulation_configure_box(byref(self), c_double(boxsize), c_int(N_root_x), c_int(N_root_y), c_int(N_root_z))
        return
   
    def configure_ghostboxes(self, N_ghost_x=0, N_ghost_y=0, N_ghost_z=0):
        """
        Initialize the ghost boxes.

        This function only needs to be called it boundary conditions other than "none" or
        "open" are used. In such a case the number of ghostboxes must be known and is set 
        with this function. 
        
        Parameters
        ----------
        N_ghost_x, N_ghost_y, N_ghost_z : int
            The number of ghost boxes in each direction. All values default to 0 (no ghost boxes).
        """
        clibrebound.N_ghost_x = c_int(N_ghost_x)
        clibrebound.N_ghost_y = c_int(N_ghost_y)
        clibrebound.N_ghost_z = c_int(N_ghost_z)
        return


# Output to file (Simulationarchive)
    def save_to_file(self, filename, interval=None, walltime=None, step=None, delete_file=False):
        """
        Saves a simulation to a file (Simulationarchive binary format). 

        If just the filename is passed, then the simulation is saved immediately. 
        You can also automate taking snapshots during a simulation. For that specify
        the interval at which outputs are required. You can do that using either
        `interval` (physical time), `walltime` (in seconds), or `step` 
        (the number of timesteps in-between snapshots). 

        
        Arguments
        ---------
        filename : str
            Filename of the binary file.
        interval : float
            Interval between outputs in code units.
        walltime : float
            Interval between outputs in wall time (seconds). 
            Useful when using IAS15 with adaptive timesteps. 
        step : int
            Interval between outputs in number of timesteps. 
            Useful when outputs need to be spaced exactly.
        delete_file: bool
            False (default) appends to file if it exists.
            True deletes the file first.
        
        Examples
        --------
        The following example shows how to output a simulation and
        then read it back in. 

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.save_to_file("test.bin")
        >>> sim2 = rebound.Simulation("test.bin")

        The following example creates a simulation, then 
        initializes the Simulationarchive and integrates
        it forward in time. 

        >>> sim = rebound.Simulation()
        >>> sim.add(m=1.)
        >>> sim.add(m=1.e-3,x=1.,vy=1.)
        >>> sim.save_to_file("sa.bin",interval=1000.)
        >>> sim.integrate(1e8)

        The Simulationarchive can later be read in using the following syntax:

        >>> sa = rebound.Simulationarchive("sa.bin")
        >>> sim = sa[0]   # get the first snapshot in the SA file (initial conditions)
        >>> sim = sa[-1]  # get the last snapshot in the SA file

        """
        if delete_file and os.path.isfile(filename):
            os.remove(filename)

        modes = sum(1 for i in [interval, walltime,step] if i is not None)

        if modes == 0:
            # Immediately save to file.
            clibrebound.reb_simulation_save_to_file(byref(self), c_char_p(filename.encode("ascii")))
        elif modes == 1:
            if delete_file:
                # reset intervals so that automate functions in C set sim->next consistently
                self.simulationarchive_auto_interval = 0
                self.simulationarchive_auto_walltime = 0
                self.simulationarchive_auto_step = 0
            if interval:
                clibrebound.reb_simulation_save_to_file_interval(byref(self), c_char_p(filename.encode("ascii")), c_double(interval))
            if walltime:
                clibrebound.reb_simulation_save_to_file_walltime(byref(self), c_char_p(filename.encode("ascii")), c_double(walltime))
            if step:
                clibrebound.reb_simulation_save_to_file_step(byref(self), c_char_p(filename.encode("ascii")), c_uint64(step))
            self.process_messages()
        else:
            raise AttributeError("Cannot specify more than one of interval, walltime, or step.")


# Integration
    def step(self):
        """
        Perform exactly one integration step with REBOUND. This function is rarely needed.
        Instead, use integrate().
        """
        clibrebound.reb_simulation_step(byref(self))
        self.process_messages()
    
    def steps(self, N_steps):
        """
        Perform exactly N_steps integration steps with REBOUND. This function is rarely needed.
        Instead, use integrate().
        """
        clibrebound.reb_simulation_steps(byref(self),c_uint(N_steps))
        self.process_messages()

    def integrate(self, tmax, exact_finish_time=1):
        """
        Main integration function. Call this function when you have setup your simulation and want to integrate it forward (or backward) in time. The function might be called many times to integrate the simulation in steps and create outputs in-between steps.
        
        Parameters
        ----------
        tmax : float
            The final time of your simulation. If the current time is 100, and tmax=200, then after the calling the integrate routine, the time has advanced to t=200. If tmax is larger than or equal to the current time, no integration will be performed.
        exact_finish_time: int, optional
            This argument determines whether REBOUND should try to finish at the exact time (tmax) you give it or if it is allowed to overshoot. Overshooting could happen if one starts at t=0, has a timestep of dt=10 and wants to integrate to tmax=25. With ``exact_finish_time=1``, the integrator will choose the last timestep such that t is exactly 25 after the integration, otherwise t=30. Note that changing the timestep does affect the accuracy of symplectic integrators negatively.
        
        Exceptions
        ----------
        Exceptions are thrown when no more particles are left in the simulation or when a generic integration error occured. 
        If you specified exit_min_distance or exit_max_distance, then additional exceptions might thrown for escaping particles or particles that undergo a clos encounter.
        
        Examples
        -------- 
        The typical usage is as follows. Note the use of ``np.linspace`` to create equally spaced outputs. 
        Using ``np.logspace`` can be used to easily produce logarithmically spaced outputs.

        >>> import numpy as np
        >>> for time in np.linspace(0,100.,10): 
        >>>     sim.integrate(time)
        >>>     perform_output(sim)
        
        """
        self.exact_finish_time = c_int(exact_finish_time)
        ret_value = clibrebound.reb_simulation_integrate(byref(self), c_double(tmax))
        if ret_value == 1:
            self.process_messages()
            raise GenericError("An error occured during the integration.")
        if ret_value == 2:
            if self._N_odes>0:
                raise NoParticles("No particles found. Will exit. Use BS integrator to integrate user-defined ODEs without any particles present.");
            else:
                raise NoParticles("No more particles left in simulation.")
        if ret_value == 3:
            raise Encounter("Two particles had a close encounter (d<exit_min_distance).")
        if ret_value == 4:
            raise Escape("A particle escaped (r>exit_max_distance).")
        if ret_value == 5:
            pass # User caused exit. Do not raise error message
        if ret_value == 6:
            raise KeyboardInterrupt
        if ret_value == 7:
            raise Collision("Two particles collided (d < r1+r2)")
        self.process_messages()

    def stop(self):
        """
        Call this function to stop an integration, for example
        from the heartbeat function.
        """
        clibrebound.reb_simulation_stop(byref(self))

    def integrator_reset(self):
        """
        Call this function to reset temporary integrator variables
        """
        clibrebound.reb_simulation_reset_integrator(byref(self))

    def integrator_synchronize(self):
        """
        Call this function if safe-mode is disabled and you need to synchronize particle positions and velocities between timesteps.
        """
        clibrebound.reb_simulation_synchronize(byref(self))
    
    def tree_update(self):
        """
        Call this function to update the tree structure manually after removing particles.
        """
        clibrebound.reb_simulation_update_tree(byref(self))
    

class timeval(Structure):
    _fields_ = [("tv_sec",c_int64),("tv_usec",c_int64)]

from .particle import Particle
from .particles import Particles


from .integrators.bs import ODE, IntegratorBS
from .integrators.whfast import IntegratorWHFast
from .integrators.whfast512 import IntegratorWHFast512
from .integrators.janus import IntegratorJanus
from .integrators.sei import IntegratorSEI
from .integrators.eos import IntegratorEOS
from .integrators.ias15 import IntegratorIAS15
from .integrators.saba import IntegratorSABA
from .integrators.mercurius import IntegratorMercurius

from .variation import Variation

# Setting up fields after class definition (because of self-reference)

class ServerData(Structure):
    _fields_ = [
            ("r", POINTER(Simulation)),
            ("r_copy", POINTER(Simulation)),
            ("port", c_int),
            ("need_copy", c_int),
            ("ready", c_int),
            # other fields not needed.
            ]

Simulation._fields_ = [
                ("t", c_double),
                ("G", c_double),
                ("softening", c_double),
                ("dt", c_double),
                ("dt_last_done", c_double),
                ("steps_done", c_uint64),
                ("N", c_uint),
                ("N_var", c_int),
                ("N_var_config", c_int),
                ("var_config", POINTER(Variation)),
                ("_var_rescale_warning", c_int),
                ("N_active", c_int),
                ("testparticle_type", c_int),
                ("testparticle_hidewarnings", c_int),
                ("_particle_lookup_table", POINTER(HashPointerPair)),
                ("hash_ctr", c_int),
                ("N_lookup", c_int),
                ("N_allocated_lookup", c_int),
                ("N_allocated", c_uint),
                ("_particles", POINTER(Particle)),
                ("gravity_cs", POINTER(Vec3dBasic)),
                ("N_allocated_gravity_cs", c_int),
                ("_tree_root", c_void_p),
                ("_tree_needs_update", c_int),
                ("opening_angle2", c_double),
                ("_status", c_int),
                ("exact_finish_time", c_int),
                ("force_is_velocity_dependent", c_uint),
                ("gravity_ignore", c_uint),
                ("_output_timing_last", c_double),
                ("save_messages", c_int),
                ("messages", c_void_p),
                ("exit_max_distance", c_double),
                ("exit_min_distance", c_double),
                ("usleep", c_double),
                ("_display_view", c_void_p),
                ("_display_data", c_void_p), # not needed from python
                ("_server_data", POINTER(ServerData)),
                ("track_energy_offset", c_int),
                ("energy_offset", c_double),
                ("walltime", c_double),
                ("walltime_last_step", c_double),
                ("walltime_last_steps", c_double),
                ("_walltime_last_steps_sum", c_double),
                ("_walltime_last_steps_N", c_int),
                ("python_unit_t",c_uint32),
                ("python_unit_l",c_uint32),
                ("python_unit_m",c_uint32),
                ("boxsize", Vec3dBasic),
                ("boxsize_max", c_double),
                ("root_size", c_double),
                ("N_root", c_int),
                ("N_root_x", c_int),
                ("N_root_y", c_int),
                ("N_root_z", c_int),
                ("N_ghost_x", c_int),
                ("N_ghost_y", c_int),
                ("N_ghost_z", c_int),
                ("collision_resolve_keep_sorted", c_int),
                ("collisions", c_void_p),
                ("N_allocated_collisions", c_int),
                ("minimum_collision_velocity", c_double),
                ("collisions_plog", c_double),
                ("max_radius", c_double*2),
                ("collisions_log_n", c_int64),
                ("_calculate_megno", c_int),
                ("_megno_Ys", c_double),
                ("_megno_Yss", c_double),
                ("_megno_cov_Yt", c_double),
                ("_megno_var_t", c_double),
                ("_megno_mean_t", c_double),
                ("_megno_mean_Y", c_double),
                ("_megno_n", c_int64),
                ("rand_seed",c_uint),
                ("simulationarchive_version", c_int),
                ("simulationarchive_auto_interval", c_double),
                ("simulationarchive_auto_walltime", c_double),
                ("simulationarchive_auto_step", c_uint64),
                ("simulationarchive_next", c_double),
                ("simulationarchive_next_step", c_uint64),
                ("_simulationarchive_filename", c_char_p),
                ("_collision", c_int),
                ("_integrator", c_int),
                ("_boundary", c_int),
                ("_gravity", c_int),
                ("ri_sei", IntegratorSEI), 
                ("ri_whfast", IntegratorWHFast),
                ("ri_whfast512", IntegratorWHFast512),
                ("ri_saba", IntegratorSABA),
                ("ri_ias15", IntegratorIAS15),
                ("ri_mercurius", IntegratorMercurius),
                ("ri_janus", IntegratorJanus),
                ("ri_eos", IntegratorEOS),
                ("ri_bs", IntegratorBS),
                ("_odes", POINTER(POINTER(ODE))),
                ("_N_odes", c_int),
                ("_N_allocated_odes", c_int),
                ("_odes_warnings", c_int),
                ("_additional_forces", CFUNCTYPE(None,POINTER(Simulation))),
                ("_pre_timestep_modifications", CFUNCTYPE(None,POINTER(Simulation))),
                ("_post_timestep_modifications", CFUNCTYPE(None,POINTER(Simulation))),
                ("_heartbeat", CFUNCTYPE(None,POINTER(Simulation))),
                ("_key_callback", CFUNCTYPE(c_int,POINTER(Simulation), c_int)),
                ("_coefficient_of_restitution", CFUNCTYPE(c_double,POINTER(Simulation), c_double)),
                ("_collision_resolve", CFUNCTYPE(c_int,POINTER(Simulation), CollisionS)),
                ("_free_particle_ap", CFUNCTYPE(None, POINTER(Particle))),
                ("_extras_cleanup", CFUNCTYPE(None, POINTER(Simulation))),
                ("extras", c_void_p),
                 ]


AFF = CFUNCTYPE(None,POINTER(Simulation))
CORFF = CFUNCTYPE(c_double,POINTER(Simulation), c_double)
COLRFF = CFUNCTYPE(c_int, POINTER(Simulation), CollisionS)
FPA = CFUNCTYPE(None, POINTER(Particle))


# Import at the end to avoid circular dependence
from . import horizons
from . import data
from .simulationarchive import Simulationarchive
