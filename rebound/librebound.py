from ctypes import *
import math
import os
import ctypes.util
import pkg_resources
import types

class ReboundModule(types.ModuleType):
    _pymodulespath = os.path.dirname(__file__)
#Find the rebound C library
    clibrebound = None
    try:
        clibrebound = CDLL(_pymodulespath+"/../librebound.so", RTLD_GLOBAL)
    except:
        print("Cannot find library 'librebound.so'.")
        raise

    AFF = CFUNCTYPE(None)
    fp = None

    @property
    def build_str():
        return c_char_p.in_dll(self.clibrebound, "build_str").value

    def status():
        """ Returns a string with a summary of the current status 
            of the simulation
            """
        s= ""
        N = get_N()
        s += "---------------------------------\n"
        s += "Rebound version:     \t" + pkg_resources.require("rebound")[0].version +"\n"
        s += "Build on:            \t" +get_build_str() + "\n"
        s += "Number of particles: \t%d\n" %N       
        s += "Simulation time:     \t%f\n" %get_t()
        if N>0:
            s += "---------------------------------\n"
            p = self.particles()
            for i in range(N):
                s += str(p[i]) + "\n"
        s += "---------------------------------"
        print(s)

# Set function pointer for additional forces
    def set_additional_forces(func):
        if(isinstance(func,types.FunctionType)):
            # Python function pointer
            global fp  # keep references
            fp = AFF(func)
            self.clibrebound.set_additional_forces(fp)
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
    def min_dt(self,t):
        return c_double.in_dll(self.clibrebound, "integrator_ias15_min_dt")

    @min_dt.setter
    def min_dt(self, value):
        c_double.in_dll(self.clibrebound, "integrator_ias15_min_dt").value = value

    @property
    def N(self):
        return c_int.in_dll(self.clibrebound,"N").value 

    @property
    def iter(self):
        return c_int.in_dll(self.clibrebound,"iter").value 

    @property
    def timing(self):
        return c_double.in_dll(self.clibrebound,"timing").value 
    
    
# MEGNO
    def init_megno(self,delta):
        self.clibrebound.tools_megno_init(c_double(delta))
    
    @property
    def megno(self):
        self.clibrebound.tools_megno.restype = c_double
        return self.clibrebound.tools_megno()
    
    @property
    def lyapunov(self):
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
        """
        if particle is not None:
            if isinstance(particle,Particle):
                self.clibrebound.particles_add(particle)
            elif isinstance(particle,str):
                self.add(horizons.getParticle(particle,**kwargs))
        else: 
            self.add(Particle(**kwargs))

# Particle getter functions
    def get_particle(self,i):
        N = get_N() 
        if i>=N:
            return None
        getp = self.clibrebound.particle_get
        getp.restype = Particle
        _p = getp(c_int(i))
        return _p

    @property
    def particles(self):
        """Return an array that points to the particle structure.
        This is a pointer and thus the contents of the array update 
        as the simulation progresses. Note that the pointer could change,
        for example when a particle is added or removed from the simulation. 
        """
        N = c_int.in_dll(self.clibrebound,"N").value 
        getp = self.clibrebound.particles_get
        getp.restype = POINTER(Particle)
        return getp()

    def remove_all_particles(self):
        self.clibrebound.particles_remove_all()


# Orbit getter
    def get_orbits(self,heliocentric=False):
        """ Returns an array of Orbits of length N-1.

        Parameters
        __________

        By default this functions returns the orbits in Jacobi coordinates. 
        Set the parameter heliocentric to True to return orbits in heliocentric coordinates.
        """
        _particles = self.particles()
        orbits = []
        for i in range(1,get_N()):
            if heliocentric:
                com = _particles[0]
            else:
                com = get_com(i)
            orbits.append(_particles[i].get_orbit(primary=com))
        return orbits

# Tools
    def move_to_com(self):
        self.clibrebound.tools_move_to_center_of_momentum()

    def reset(self):
        debug.reset_debug()
        self.clibrebound.reset()

    def set_integrator_whfast_corrector(self,on=11):
        c_int.in_dll(self.clibrebound, "integrator_whfast_corrector").value = on

    def set_integrator(self,integrator="IAS15"):
        if isinstance(integrator, int):
            self.clibrebound.integrator_set(c_int(integrator))
            return
        if isinstance(integrator, basestring):
            debug.integrator_fullname = integrator
            debug.integrator_package = "REBOUND"
            if integrator.lower() == "ias15":
                set_integrator(0)
                return
            if integrator.lower() == "whfast":
                set_integrator(1)
                set_integrator_whfast_corrector(11)
                return
            if integrator.lower() == "whfast-nocor":
                set_integrator(1)
                set_integrator_whfast_corrector(0)
                return
            if integrator.lower() == "sei":
                set_integrator(2)
                return
            if integrator.lower() == "wh":
                set_integrator(3)
                return
            if integrator.lower() == "leapfrog":
                set_integrator(4)
                return
            if integrator.lower() == "leap-frog":
                set_integrator(4)
                return
            if integrator.lower() == "hybrid":
                set_integrator(5)
                return
            if integrator.lower() == "mercury":
                debug.integrator_package = "MERCURY"
                return
            if integrator.lower() == "swifter-whm":
                debug.integrator_package = "SWIFTER"
                return
            if integrator.lower() == "swifter-symba":
                debug.integrator_package = "SWIFTER"
                return
            if integrator.lower() == "swifter-helio":
                debug.integrator_package = "SWIFTER"
                return
            if integrator.lower() == "swifter-tu4":
                debug.integrator_package = "SWIFTER"
                return
        raise ValueError("Warning. Intergrator not found.")

    def set_integrator_whfast_persistent_particles(self,is_per=0):
        if isinstance(is_per, int):
            c_int.in_dll(self.clibrebound, "integrator_whfast_persistent_particles").value = is_per
            return
        raise ValueError("Expecting integer.")

    def set_integrator_whfast_synchronize_manually(self,synchronize_manually=0):
        if isinstance(synchronize_manually, int):
            c_int.in_dll(self.clibrebound, "integrator_whfast_synchronize_manually").value = synchronize_manually
            return
        raise ValueError("Expecting integer.")

    def set_force_is_velocitydependent(self,force_is_velocitydependent=1):
        if isinstance(force_is_velocitydependent, int):
            c_int.in_dll(self.clibrebound, "integrator_force_is_velocitydependent").value = force_is_velocitydependent
            return
        raise ValueError("Expecting integer.")

# Input/Output routines

    def output_binary(self,filename):
        self.clibrebound.output_binary(c_char_p(filename))
    save=output_binary
        
    def input_binary(self,filename):
        self.clibrebound.input_binary(c_char_p(filename))
    load=input_binary
        

# Integration
    def step(self):
        self.clibrebound.rebound_step()

    def integrate(self,tmax,exactFinishTime=1,keepSynchronized=0,maxR=0.,minD=0.):
        if debug.integrator_package =="REBOUND":
            self.clibrebound.integrate.restype = c_int
            ret_value = self.clibrebound.integrate(c_double(tmax),c_int(exactFinishTime),c_int(keepSynchronized),c_double(maxR),c_double(minD))
            if ret_value == 1:
                raise NoParticleLeft("No more particles left in simulation.")
            if ret_value == 2:
                raise ParticleEscaping("At least one particle has a radius > maxR.")
            if ret_value == 3:
                raise CloseEncounter(c_int.in_dll(self.clibrebound, "closeEncounterPi").value,
                                     c_int.in_dll(self.clibrebound, "closeEncounterPj").value)
        else:
            debug.integrate_other_package(tmax,exactFinishTime,keepSynchronized)

    @property
    def com(self,last=None):
        """Returns the center of momentum for all particles in the simulation"""
        m = 0.
        x = 0.
        y = 0.
        z = 0.
        vx = 0.
        vy = 0.
        vz = 0.
        ps = self.particles()    # particle pointer
        if last is not None:
            last = min(last,get_N())
        else:
            last = get_N()
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

# Helper functions
TWOPI = 2.*math.pi
def mod2pi(f):
    """Returns the angle f modulo 2 pi."""
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def eccentricAnomaly(e,M):
    """Returns the eccentric anomaly given the eccentricity and mean anomaly of a Keplerian orbit.

    Keyword arguments:
    e -- the eccentricity
    M -- the mean anomaly
    """
    if e<1.:
        E = M if e<0.8 else math.pi
        
        F = E - e*math.sin(M) - M
        for i in range(100):
            E = E - F/(1.0-e*math.cos(E))
            F = E - e*math.sin(E) - M
            if math.fabs(F)<1e-16:
                break
        E = mod2pi(E)
        return E
    else:
        E = M 
        
        F = E - e*math.sinh(E) - M
        for i in range(100):
            E = E - F/(1.0-e*math.cosh(E))
            F = E - e*math.sinh(E) - M
            if math.fabs(F)<1e-16:
                break
        E = mod2pi(E)
        return E



# Import at the end to avoid circular dependence
from .particle import *
from . import horizons
from . import debug
