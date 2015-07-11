from ctypes import *
import math
import os
import ctypes.util
try:
    import pkg_resources
except: 
    # Fails on python3, but not important
    pass
import types
        
INTEGRATORS = {"ias15": 0, "whfast": 1, "sei": 2, "wh": 3, "leapfrog": 4, "hybrid": 5, "none": 6}

class ReboundModule(types.ModuleType):
    _pymodulespath = os.path.dirname(__file__)
    clibrebound = None   # Class variable
    try:
        clibrebound = CDLL(_pymodulespath+"/../librebound.so", RTLD_GLOBAL)
    except:
        print("Cannot find library 'librebound.so'.")
        raise

    AFF = CFUNCTYPE(None)
    afp = None # additional forces pointer
    ptmp = None # post timestep modifications pointer 
    _units = {'length':None, 'time':None, 'mass':None}

# Status functions
    @property
    def build_str(self):
        return str(c_char_p.in_dll(self.clibrebound, "build_str").value)

    def status(self):
        """ Returns a string with a summary of the current status 
            of the simulation
            """
        s= ""
        s += "---------------------------------\n"
        try:
            s += "Rebound version:     \t" + pkg_resources.require("rebound")[0].version +"\n"
        except:
            # Fails on python3, but not important
            pass
        s += "Build on:            \t" + self.build_str + "\n"
        s += "Number of particles: \t%d\n" %self.N       
        s += "Simulation time:     \t%f\n" %self.t
        if self.N>0:
            s += "---------------------------------\n"
            for p in self.particles:
                s += str(p) + "\n"
        s += "---------------------------------"
        print(s)

# Set function pointer for additional forces
    @property
    def additional_forces(self):
        return self.afp   # might not be needed

    @additional_forces.setter
    def additional_forces(self, func):
        if(isinstance(func,types.FunctionType)):
            # Python function pointer
            self.afp = self.AFF(func)
            self.clibrebound.set_additional_forces(self.afp)
        else:
            # C function pointer
            self.clibrebound.set_additional_forces_with_parameters(func)
            self.afp = "C function pointer value currently not accessible from python.  Edit librebound.py"

    @property
    def post_timestep_modifications(self):
        return self.ptmp

    @post_timestep_modifications.setter
    def post_timestep_modifications(self, func):
        if(isinstance(func, types.FunctionType)):
            # Python function pointer
            self.ptmp = self.AFF(func)
            self.clibrebound.set_post_timestep_modifications(self.ptmp)
        else:
            # C function pointer
            self.clibrebound.set_post_timestep_modifications_with_parameters(func)
            self.ptmp = "C function pointer value currently not accessible from python.  Edit librebound.py" 

# Setter/getter of parameters and constants
    @property
    def G(self):
        return c_double.in_dll(self.clibrebound, "G").value

    @G.setter
    def G(self, value):
        c_double.in_dll(self.clibrebound, "G").value = value

    @property
    def dt(self):
        return c_double.in_dll(self.clibrebound, "dt").value

    @dt.setter
    def dt(self, value):
        c_double.in_dll(self.clibrebound, "dt").value = value

    @property
    def t(self):
        return c_double.in_dll(self.clibrebound, "t").value

    @t.setter
    def t(self, value):
        c_double.in_dll(self.clibrebound, "t").value = value

    @property
    def min_dt(self):
        return c_double.in_dll(self.clibrebound, "integrator_ias15_min_dt").value

    @min_dt.setter
    def min_dt(self, value):
        c_double.in_dll(self.clibrebound, "integrator_ias15_min_dt").value = value
    
    @property
    def integrator_ias15_epsilon(self):
        return c_double.in_dll(self.clibrebound, "integrator_ias15_epsilon").value

    @integrator_ias15_epsilon.setter
    def integrator_ias15_epsilon(self, value):
        c_double.in_dll(self.clibrebound, "integrator_ias15_epsilon").value = value


    @property
    def N(self):
        return c_int.in_dll(self.clibrebound, "N").value 
   
    @N.setter
    def N(self, value):
        c_int.in_dll(self.clibrebound, "N").value = value

    @property
    def N_active(self):
        return c_int.in_dll(self.clibrebound, "N_active").value

    @N_active.setter
    def N_active(self, value):
        c_int.in_dll(self.clibrebound, "N_active").value = value

    @property
    def integrator(self):
        i = c_int.in_dll(self.clibrebound, "integrator").value
        for name, _i in INTEGRATORS.items():
            if i==_i:
                return name
        return i

    @integrator.setter
    def integrator(self, value):
        if isinstance(value, int):
            self.clibrebound.set_integrator(c_int(value))
        elif isinstance(value, basestring):
            debug.integrator_fullname = value
            debug.integrator_package = "REBOUND"
            value = value.lower()
            if value in INTEGRATORS: 
                self.integrator = INTEGRATORS[value]
            elif value.lower() == "mercury":
                debug.integrator_package = "MERCURY"
            elif value.lower() == "swifter-whm":
                debug.integrator_package = "SWIFTER"
            elif value.lower() == "swifter-symba":
                debug.integrator_package = "SWIFTER"
            elif value.lower() == "swifter-helio":
                debug.integrator_package = "SWIFTER"
            elif value.lower() == "swifter-tu4":
                debug.integrator_package = "SWIFTER"
            else:
                raise ValueError("Warning. Integrator not found.")

    @property
    def force_is_velocitydependent(self):
        return c_int.in_dll(self.clibrebound, "integrator_force_is_velocitydependent").value

    @force_is_velocitydependent.setter
    def force_is_velocitydependent(self, value):
        if isinstance(value, int):
            c_int.in_dll(self.clibrebound, "integrator_force_is_velocitydependent").value = value
            return
        raise ValueError("Expecting integer.")
    
# Units

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, newunits):
        newunits = self.check_units(newunits)        
        if self.particles: # some particles are loaded
            raise Exception("Error:  You cannot set the units after populating the particles array.  See Units.ipynb in python_tutorials.")
        self.update_units(newunits) 

    def check_units(self, newunits):   
        if len(newunits) is not 3:
            raise Exception("Error: Need to pass exactly 3 units for length, time, and mass (any order), see python_tutorials/Units.ipynb")
        
        l_unit = t_unit = m_unit = None
        for unit in newunits:
            unit = unit.lower()
            if unit in lengths_SI:
                l_unit = unit
            if unit in times_SI:
                t_unit = unit
            if unit in masses_SI:
                m_unit = unit

        if l_unit is None or t_unit is None or m_unit is None:
            raise Exception("Error: Need to assign rebound.units a tuple consisting of 3 units for length, time, and mass (any order).  See python_tutorials/Units.ipynb.  If you passed such a tuple, at least one of your units isn't in our list.  Please update the dictionaries at the bottom of rebound/rebound/librebound.py and send a pull request!")

        return (l_unit, t_unit, m_unit)

    def update_units(self, newunits): 
        self._units['length'] = newunits[0] 
        self._units['time'] = newunits[1] 
        self._units['mass'] = newunits[2] 
        self.G = self.convert_G(self._units['length'], self._units['time'], self._units['mass'])

    def convert_particle_units(self, *args): 
        new_l, new_t, new_m = self.check_units(args)
        for p in self.particles:
            self.convert_p(p, self._units['length'], self._units['time'], self._units['mass'], new_l, new_t, new_m)
        self.update_units((new_l, new_t, new_m))

    def convert_p(self, p, old_l, old_t, old_m, new_l, new_t, new_m):
        p.m = self.convert_mass(p.m, old_m, new_m)
        p.x = self.convert_length(p.x, old_l, new_l) 
        p.y = self.convert_length(p.y, old_l, new_l)
        p.z = self.convert_length(p.z, old_l, new_l)
        p.vx = self.convert_vel(p.vx, old_l, old_t, new_l, new_t)
        p.vy = self.convert_vel(p.vy, old_l, old_t, new_l, new_t)
        p.vz = self.convert_vel(p.vz, old_l, old_t, new_l, new_t)
        p.ax = self.convert_acc(p.ax, old_l, old_t, new_l, new_t)
        p.ay = self.convert_acc(p.ay, old_l, old_t, new_l, new_t)
        p.az = self.convert_acc(p.az, old_l, old_t, new_l, new_t)
        return p

    def convert_mass(self, mass, old_m, new_m):
        return mass*masses_SI[old_m]/masses_SI[new_m]

    def convert_length(self, length, old_l, new_l):
        return length*lengths_SI[old_l]/lengths_SI[new_l]

    def convert_vel(self, vel, old_l, old_t, new_l, new_t):
        in_SI=vel*lengths_SI[old_l]/times_SI[old_t]
        return in_SI*times_SI[new_t]/lengths_SI[new_l]
    
    def convert_acc(self, acc, old_l, old_t, new_l, new_t):
        in_SI=acc*lengths_SI[old_l]/times_SI[old_t]**2
        return in_SI*times_SI[new_t]**2/lengths_SI[new_l]

    def convert_G(self, new_l, new_t, new_m):
        return G_SI*masses_SI[new_m]*times_SI[new_t]**2/lengths_SI[new_l]**3

       
# Debug tools
    @property
    def iter(self):
        return c_int.in_dll(self.clibrebound,"iter").value 

    @property
    def timing(self):
        return c_double.in_dll(self.clibrebound,"timing").value 
   
    
# MEGNO
    def init_megno(self, delta):
        self.clibrebound.tools_megno_init(c_double(delta))
    
    def calculate_megno(self):
        self.clibrebound.tools_megno.restype = c_double
        return self.clibrebound.tools_megno()
    
    def calculate_lyapunov(self):
        self.clibrebound.tools_lyapunov.restype = c_double
        return self.clibrebound.tools_lyapunov()
    
    @property
    def N_megno(self):
        return c_int.in_dll(self.clibrebound,"N_megno").value 
    
# Particle add function, used to be called particle_add() and add_particle() 
    def add(self, particle=None, **kwargs):   
        """Adds a particle to REBOUND. Accepts one of the following:
        1) A single Particle structure.
        2) The particle's mass and a set of cartesian coordinates: m,x,y,z,vx,vy,vz.
        3) The primary as a Particle structure, the particle's mass and a set of orbital elements primary,a,anom,e,omega,inv,Omega,MEAN (see kepler_particle() for the definition of orbital elements). 
        4) A name of an object (uses NASA Horizons to look up coordinates)
        5) A list of particles or names.
        """
        if particle is not None:
            if isinstance(particle, Particle):
                if kwargs == {}: # copy particle
                    self.clibrebound.particles_add(particle)
                else: # use particle as primary
                    self.add(Particle(primary=particle, **kwargs))
            elif isinstance(particle, list):
                for p in particle:
                    self.add(p)
            elif isinstance(particle,str):
                if None in self.units.values():
                    self.units = ('AU', 'yr2pi', 'Msun')
                self.add(horizons.getParticle(particle,**kwargs))
                self.convert_p(self.particles[-1], 'km', 's', 'kg', self._units['length'], self._units['time'], self._units['mass'])
        else: 
            self.add(Particle(**kwargs))

# Particle getter functions
    @property
    def particles(self):
        """Return an array that points to the particle structure.
        This is an array of pointers and thus the contents of the array update 
        as the simulation progresses. Note that the pointers could change,
        for example when a particle is added or removed from the simulation. 
        """
        ps = []
        N = c_int.in_dll(self.clibrebound,"N").value 
        getp = self.clibrebound.get_particles
        getp.restype = POINTER(Particle)
        ps_a = getp()
        for i in range(0,N):
            ps.append(ps_a[i])
        return ps

    @particles.deleter
    def particles(self):
        self.clibrebound.particles_remove_all()

    def remove_particle(self, index=None, ID=None, keepSorted=1):
        """ Removes a particle from the simulation.

        Parameters

        ----------

        Either the index in the particles array to remove, or the ID of the particle to
        remove.  The keepSorted flag ensures the particles array remains sorted
        in order of increasing IDs.  One might set it to zero in cases with many
        particles and many removals to speed things up.
        """
        if index is not None:
            success = self.clibrebound.particles_remove(c_int(index), keepSorted)
            if not success:
                print("Index %d passed to remove_particle was out of range (N=%d). Did not remove particle.\n"%(index, self.N))
            return
        if ID is not None:
            success = self.clibrebound.particles_remove_ID(c_int(ID), keepSorted)
            if not success:
                print("ID %d passed to remove_particle was not found.  Did not remove particle.\n"%(ID))


# Orbit calculation
    def calculate_orbits(self, heliocentric=False):
        """ Returns an array of Orbits of length N-1.

        Parameters
        __________

        By default this functions returns the orbits in Jacobi coordinates. 
        Set the parameter heliocentric to True to return orbits in heliocentric coordinates.
        """
        _particles = self.particles
        orbits = []
        for i in range(1,self.N):
            if heliocentric:
                com = _particles[0]
            else:
                com = self.calculate_com(i)
            orbits.append(_particles[i].calculate_orbit(primary=com))
        return orbits

# COM calculation 
    def calculate_com(self, last=None):
        """Returns the center of momentum for all particles in the simulation"""
        m = 0.
        x = 0.
        y = 0.
        z = 0.
        vx = 0.
        vy = 0.
        vz = 0.
        ps = self.particles    # particle pointer
        if last is not None:
            last = min(last, self.N)
        else:
            last = self.N
        for i in range(last):
            m  += ps[i].m
            x  += ps[i].x*ps[i].m
            y  += ps[i].y*ps[i].m
            z  += ps[i].z*ps[i].m
            vx += ps[i].vx*ps[i].m
            vy += ps[i].vy*ps[i].m
            vz += ps[i].vz*ps[i].m
        if m>0.:
            x /= m
            y /= m
            z /= m
            vx /= m
            vy /= m
            vz /= m
        return Particle(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    

# Tools
    def move_to_com(self):
        self.clibrebound.tools_move_to_center_of_momentum()
    
    def calculate_energy(self):
        self.clibrebound.tools_energy.restype = c_double
        return self.clibrebound.tools_energy()

    def reset(self):
        self._units = {'length':None, 'time':None, 'mass':None}
        debug.reset_debug()
        self.clibrebound.reset()


# Input/Output routines
    def save(self, filename):
        self.clibrebound.output_binary(c_char_p(filename.encode("ascii")))
        
    def load(self, filename):
        self.clibrebound.input_binary(c_char_p(filename.encode("ascii")))
        
# Integrator Flags
    @property 
    def integrator_whfast_corrector(self):
        return c_int.in_dll(self.clibrebound, "integrator_whfast_corrector").value
    
    @integrator_whfast_corrector.setter 
    def integrator_whfast_corrector(self, value):
        c_int.in_dll(self.clibrebound, "integrator_whfast_corrector").value = value

    @property
    def integrator_whfast_safe_mode(self):
        return c_int.in_dll(self.clibrebound, "integrator_whfast_safe_mode").value

    @integrator_whfast_safe_mode.setter
    def integrator_whfast_safe_mode(self, value):
        c_int.in_dll(self.clibrebound, "integrator_whfast_safe_mode").value = value

    @property
    def recalculate_jacobi_this_timestep(self):
        return c_int.in_dll(self.clibrebound, "integrator_whfast_recalculate_jacobi_this_timestep").value

    @recalculate_jacobi_this_timestep.setter
    def recalculate_jacobi_this_timestep(self, value):
        c_int.in_dll(self.clibrebound, "integrator_whfast_recalculate_jacobi_this_timestep").value = value
    
# Integration

    def step(self, do_timing = 1):
        self.clibrebound.rebound_step(c_int(do_timing))

    def integrate(self, tmax, exact_finish_time=1, maxR=0., minD=0.):
        if debug.integrator_package =="REBOUND":
            self.clibrebound.rebound_integrate.restype = c_int
            ret_value = self.clibrebound.rebound_integrate(c_double(tmax),c_int(exact_finish_time),c_double(maxR),c_double(minD))
            if ret_value == 1:
                raise self.NoParticleLeft("No more particles left in simulation.")
            if ret_value == 2:
                raise self.ParticleEscaping(c_int.in_dll(self.clibrebound, "escapedParticle").value, self.t)
            if ret_value == 3:
                raise self.CloseEncounter(c_int.in_dll(self.clibrebound, "closeEncounterPi").value,
                                     c_int.in_dll(self.clibrebound, "closeEncounterPj").value, self.t)
        else:
            debug.integrate_other_package(tmax,exact_finish_time)

    def integrator_synchronize(self):
        self.clibrebound.integrator_synchronize()

# Exceptions    
    class CloseEncounter(Exception):
        def __init__(self, index1, index2, t):
                self.index1 = index1
                self.index2 = index2
                self.t = t
        def __str__(self):
                return "A close encounter occured at time %f between particles with indices %d and %d."%(self.t, self.index1,self.index2)

    class ParticleEscaping(Exception):
        def __init__(self, index, t):
            self.index = index
            self.t = t
        def __str__(self):
            return "At least one particle's distance > maxR at time %f (index = %d)."%(self.t, self.index)

    class NoParticleLeft(Exception):
        pass

# Unit constants
G_SI = 6.674e-11
times_SI = {'s':1.,
    'hr':3600.,
    'yr':31557600., # Julian year (exact)
    'jyr':31557600.,
    'sidereal_yr':31558149.7635,
    'yr2pi':math.sqrt(149597870700.**3/1.3271244004193938e20), # chosen to make G=1
    'kyr':31557600.*1.e3,
    'myr':31557600.*1.e6,
    'gyr':31557600.*1.e9}
lengths_SI =  {'m':1.,
    'cm':0.01,
    'km':1000.,
    'au':149597870700.}

    #What we measure accurately is GM, so set mass units such that G*M gives the value of GM in horizons.py (in the list at the end of horizons.py, the NAIF codes ending in 99 refer to the planets, single digits to the total mass of the planet plus its moons).  Have to multiply by 10**9 since that list has G in kg^-1km^3/s^2 and we use SI.

masses_SI = {'kg':1.,
    'msun':1.3271244004193938E+11/G_SI*10**9,
    'mmercury':2.2031780000000021E+04/G_SI*10**9,
    'mvenus':3.2485859200000006E+05/G_SI*10**9,
    'mearth':3.9860043543609598E+05/G_SI*10**9,
    'mmars':4.282837362069909E+04/G_SI*10**9,
    'mjupiter':1.266865349218008E+08/G_SI*10**9,
    'msaturn':3.793120749865224E+07/G_SI*10**9,
    'muranus':5.793951322279009E+06/G_SI*10**9,
    'mneptune':6.835099502439672E+06/G_SI*10**9,
    'mpluto':8.696138177608748E+02/G_SI*10**9}



# Import at the end to avoid circular dependence
from .particle import *
from . import horizons
from . import debug
