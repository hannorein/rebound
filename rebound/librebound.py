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
    fp = None
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
        return self.fp   # might not be needed

    @additional_forces.setter
    def additional_forces(self,func):
        if(isinstance(func,types.FunctionType)):
            # Python function pointer
            self.fp = self.AFF(func)
            self.clibrebound.set_additional_forces(self.fp)
        else:
            # C function pointer
            self.clibrebound.set_additional_forces_with_parameters(func)

# Setter/getter of parameters and constants
    @property
    def G(self):
        return c_double.in_dll(self.clibrebound, "G").value

    @G.setter
    def G(self,value):
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
        return c_int.in_dll(self.clibrebound,"N").value 
    
    @property
    def N_active(self):
        return c_int.in_dll(self.clibrebound, "N_active").value

    @N_active.setter
    def N_active(self, value):
        c_int.in_dll(self.clibrebound, "N_active").value = value

    @property 
    def integrator_whfast_corrector(self):
        return c_int.in_dll(self.clibrebound, "integrator_whfast_corrector").value
    
    @integrator_whfast_corrector.setter 
    def integrator_whfast_corrector(self, value):
        c_int.in_dll(self.clibrebound, "integrator_whfast_corrector").value = value

    @property
    def integrator(self):
        i = c_int.in_dll(self.clibrebound, "integrator").value
        for name, _i in list.iteritems():
            if i==_i:
                return name
        return i

    @integrator.setter
    def integrator(self,value):
        if isinstance(value, int):
            self.clibrebound.integrator_set(c_int(value))
        elif isinstance(value, basestring):
            debug.integrator_fullname = value
            debug.integrator_package = "REBOUND"
            value = value.lower()
            if value in INTEGRATORS: 
                self.integrator = INTEGRATORS[value]
                self.integrator_whfast_corrector = 11
            elif value.lower() == "whfast-nocor":
                self.integrator = INTEGRATORS["whfast"]
                self.integrator_whfast_corrector = 0
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
                raise ValueError("Warning. Intergrator not found.")

    @property
    def force_is_velocitydependent(self):
        return c_int.in_dll(self.clibrebound, "integrator_force_is_velocitydependent").value

    @force_is_velocitydependent.setter
    def force_is_velocitydependent(self,value):
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
        l_unit, t_unit, m_unit = newunits
        if self.particles: # some particles are loaded
            print("Error:  You must initialize the units before populating the particles array.  See Units.ipynb in python_tutorials.")
            import sys
            sys.exit()

        self._units['length'] = l_unit.lower()
        self._units['time'] = t_unit.lower()
        self._units['mass'] = m_unit.lower()
        self.G = self.convert_G()
    
    def convert_units(self, l_unit, t_unit, m_unit): 
        for p in self.particles:
            p.m = self.convert_mass(p.m, m_unit.lower())
            p.x = self.convert_length(p.x, l_unit.lower()) 
            p.y = self.convert_length(p.y, l_unit.lower()) 
            p.z = self.convert_length(p.z, l_unit.lower())
            p.vx = self.convert_vel(p.vx, l_unit.lower(), t_unit.lower()) 
            p.vy = self.convert_vel(p.vy, l_unit.lower(), t_unit.lower()) 
            p.vz = self.convert_vel(p.vz, l_unit.lower(), t_unit.lower())
            p.ax = self.convert_vel(p.ax, l_unit.lower(), t_unit.lower())
            p.ay = self.convert_vel(p.ay, l_unit.lower(), t_unit.lower())
            p.az = self.convert_vel(p.az, l_unit.lower(), t_unit.lower())

        self._units['length'] = l_unit.lower()
        self._units['time'] = t_unit.lower()
        self._units['mass'] = m_unit.lower()

        self.G = self.convert_G()

    def convert_mass(self, mass, m_unit):
        return mass*masses_SI[self._units['mass']]/masses_SI[m_unit]

    def convert_length(self, length, l_unit):
        return length*lengths_SI[self._units['length']]/lengths_SI[l_unit]

    def convert_vel(self, vel, l_unit, t_unit):
        in_SI=vel*lengths_SI[self._units['length']]/times_SI[self._units['time']]
        return in_SI*times_SI[t_unit]/lengths_SI[l_unit]

    def convert_acc(self, acc, l_unit, t_unit):
        in_SI=acc*lengths_SI[self._units['length']]/times_SI[self._units['time']]**2
        return in_SI*times_SI[t_unit]**2/lengths_SI[l_unit]
    
    def convert_G(self):
        return G_SI*masses_SI[self._units['mass']]*times_SI[self._units['time']]**2/lengths_SI[self._units['length']]**3

       
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
            if isinstance(particle,Particle):
                self.clibrebound.particles_add(particle)
            elif isinstance(particle,list):
                for p in particle:
                    self.add(p)
            elif isinstance(particle,str):
                self.add(horizons.getParticle(particle,**kwargs))
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
        getp = self.clibrebound.particles_get
        getp.restype = POINTER(Particle)
        ps_a = getp()
        for i in range(0,N):
            ps.append(ps_a[i])
        return ps

    @particles.deleter
    def particles(self):
        self.clibrebound.particles_remove_all()


# Orbit calculation
    def calculate_orbits(self,heliocentric=False):
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
    def calculate_com(self,last=None):
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
            last = min(last,self.N)
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
        debug.reset_debug()
        self.clibrebound.reset()


# Input/Output routines
    def save(self, filename):
        self.clibrebound.output_binary(c_char_p(filename.encode("ascii")))
        
    def load(self, filename):
        self.clibrebound.input_binary(c_char_p(filename.encode("ascii")))
        

# Integration
    def step(self):
        self.clibrebound.rebound_step()

    def integrate(self, tmax, exactFinishTime=1, keepSynchronized=0, particlesModified=1, maxR=0., minD=0.):
        if debug.integrator_package =="REBOUND":
            self.clibrebound.integrate.restype = c_int
            ret_value = self.clibrebound.integrate(c_double(tmax),c_int(exactFinishTime),c_int(keepSynchronized),c_int(particlesModified),c_double(maxR),c_double(minD))
            if ret_value == 1:
                raise self.NoParticleLeft("No more particles left in simulation.")
            if ret_value == 2:
                raise self.ParticleEscaping("At least one particle has a radius > maxR.")
            if ret_value == 3:
                raise self.CloseEncounter(c_int.in_dll(self.clibrebound, "closeEncounterPi").value,
                                     c_int.in_dll(self.clibrebound, "closeEncounterPj").value)
        else:
            debug.integrate_other_package(tmax,exactFinishTime,keepSynchronized)

# Exceptions    
    class CloseEncounter(Exception):
        def __init__(self, id1, id2):
                self.id1 = id1
                self.id2 = id2
        def __str__(self):
                return "A close encounter occured between particles %d and %d."%(self.id1,self.id2)

    class ParticleEscaping(Exception):
        pass

    class NoParticleLeft(Exception):
        pass

# Unit constants
G_SI = 6.674e-11
times_SI = {'s':1.,
    'hr':3600.,
    'yr':31557600., # Julian year (exact)
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
