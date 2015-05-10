from ctypes import *
import math
import os
import tempfile
import shutil
import time
import ctypes.util
import pkg_resources

_pymodulespath = os.path.dirname(__file__)
#Find the rebound C library
try:
    clibrebound = CDLL(_pymodulespath+"/../librebound.so", RTLD_GLOBAL)
except:
    try:
        clibrebound = CDLL(_pymodulespath + '/../shared/librebound/librebound.so', RTLD_GLOBAL)
    except:
        print("Cannot find library 'librebound.so'.")
        raise
    pass

#Make changes for python 2 and 3 compatibility
try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x

TINY=1.e-308
MIN_REL_ERROR = 1.e-12


def get_build_str():
    return c_char_p.in_dll(clibrebound, "build_str").value


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
        p = get_particles()
        for i in range(N):
            s += str(p[i]) + "\n"
    s += "---------------------------------"
    print(s)
# Alias
get_status = status

# Set function pointer for additional forces

AFF = CFUNCTYPE(None)
fp = None
def set_additional_forces(func):
    global fp  # keep references
    fp = AFF(func)
    clibrebound.set_additional_forces(fp)

# Setter/getter of parameters and constants
def set_G(G):
    c_double.in_dll(clibrebound, "G").value = G

def get_G():
    return c_double.in_dll(clibrebound, "G").value

def set_dt(dt):
    c_double.in_dll(clibrebound, "dt").value = dt

def get_dt():
    return c_double.in_dll(clibrebound, "dt").value

def set_t(t):
    c_double.in_dll(clibrebound, "t").value = t

def set_min_dt(t):
    c_double.in_dll(clibrebound, "integrator_ias15_min_dt").value = t

def get_t():
    return c_double.in_dll(clibrebound, "t").value

def init_megno(delta):
    clibrebound.tools_megno_init(c_double(delta))

def get_megno():
    clibrebound.tools_megno.restype = c_double
    return clibrebound.tools_megno()

def get_lyapunov():
    clibrebound.tools_lyapunov.restype = c_double
    return clibrebound.tools_lyapunov()

def get_N():
    return c_int.in_dll(clibrebound,"N").value 

def get_N_megno():
    return c_int.in_dll(clibrebound,"N_megno").value 

def get_iter():
    return c_int.in_dll(clibrebound,"iter").value 

def get_timing():
    return c_double.in_dll(clibrebound,"timing").value 

# Setter functions, particle data
def set_particles(_particles):
    c_int.in_dll(clibrebound,"N").value = len(_particles)
    arr = (Particle * len(_particles))(*_particles)
    clibrebound.setp(byref(arr))

def add_particle(particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None, date=None):   
    """Adds a particle to REBOUND. Accepts one of the following four sets of arguments:
    1) A single Particle structure.
    2) A list of Particle structures.
    3) The particle's mass and a set of cartesian coordinates: m,x,y,z,vx,vy,vz.
    3) The primary as a Particle structure, the particle's mass and a set of orbital elements primary,a,anom,e,omega,inv,Omega,MEAN (see kepler_particle() for the definition of orbital elements). 
    """
    if particle is None:
        particle = Particle(**locals())
    if isinstance(particle,list):
        for particle in particle:
            clibrebound.particles_add(particle)
    elif isinstance(particle,str):
        clibrebound.particles_add(horizons.getParticle(**locals()))
    else:
        clibrebound.particles_add(particle)

# Aliases
add_particles = add_particle
add = add_particle

# Particle getter functions
def get_particle(i):
    N = get_N() 
    if i>=N:
        return None
    getp = clibrebound.particle_get
    getp.restype = Particle
    _p = getp(c_int(i))
    return _p

def get_particles_array():
    N = get_N() 
    _particles = []
    for i in range(0,N):
        _particles.append(particle_get(i))
    return _particles


def get_particles():
    """Return an array that points to the particle structure.
    This is a pointer and thus the contents of the array update 
    as the simulation progresses. Note that the pointer could change,
    for example when a particle is added or removed from the simulation. 
    In that case, get a fresh pointer with get_particles().

    particles() is an alias of this function.
    """
    N = c_int.in_dll(clibrebound,"N").value 
    getp = clibrebound.particles_get
    getp.restype = POINTER(Particle)
    return getp()
# Alias
particles = get_particles

# Orbit getter
def get_orbits(heliocentric=False):
    """ Returns an array of Orbits of length N-1.

    Parameters
    __________

    By default this functions returns the orbits in Jacobi coordinates. 
    Set the parameter heliocentric to True to return orbits in heliocentric coordinates.
    """
    _particles = get_particles()
    orbits = []
    for i in range(1,get_N()):
        if heliocentric:
            com = _particles[0]
        else:
            com = get_com(i)
        orbits.append(_particles[i].get_orbit(primary=com))
    return orbits

# Tools
def move_to_center_of_momentum():
    clibrebound.tools_move_to_center_of_momentum()
move_to_barycentric_fram = move_to_center_of_momentum

tmpdir = None
def reset():
    global tmpdir
    if tmpdir:
        shutil.rmtree(tmpdir)
        tmpdir = None
    clibrebound.reset()

def set_integrator_whfast_corrector(on=11):
    c_int.in_dll(clibrebound, "integrator_whfast_corrector").value = on

#REMOVE Following variables
_integrator_package = "REBOUND"
_integrator_debug = ""    

def set_integrator(integrator="IAS15"):
    global _integrator_package
    global _integrator_debug
    if isinstance(integrator, int):
        clibrebound.integrator_set(c_int(integrator))
        return
    if isinstance(integrator, basestring):
        _integrator_debug = integrator
        _integrator_package = "REBOUND"
        if integrator.lower() == "ias15":
            set_integrator(0)
            return
        if integrator[0:7].lower() == "mikkola":
            set_integrator(1)
            set_integrator_whfast_corrector(0)
            if integrator[8:] == "cor3":
                set_integrator_whfast_corrector(3)
            if integrator[8:] == "cor5":
                set_integrator_whfast_corrector(5)
            if integrator[8:] == "cor7":
                set_integrator_whfast_corrector(7)
            if integrator[8:] == "cor11":
                set_integrator_whfast_corrector(11)
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
            _integrator_package = "MERCURY"
            return
        if integrator.lower() == "swifter-whm":
            _integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-symba":
            _integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-helio":
            _integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-tu4":
            _integrator_package = "SWIFTER"
            return
    raise ValueError("Warning. Intergrator not found.")

def set_integrator_whfast_persistent_particles(is_per=0):
    if isinstance(is_per, int):
        c_int.in_dll(clibrebound, "integrator_whfast_persistent_particles").value = is_per
        return
    raise ValueError("Expecting integer.")

def set_integrator_whfast_synchronize_manually(synchronize_manually=0):
    if isinstance(synchronize_manually, int):
        c_int.in_dll(clibrebound, "integrator_whfast_synchronize_manually").value = synchronize_manually
        return
    raise ValueError("Expecting integer.")

def set_force_is_velocitydependent(force_is_velocitydependent=1):
    if isinstance(force_is_velocitydependent, int):
        c_int.in_dll(clibrebound, "integrator_force_is_velocitydependent").value = force_is_velocitydependent
        return
    raise ValueError("Expecting integer.")

# Input/Output routines

def output_binary(filename):
    clibrebound.output_binary(c_char_p(filename))
    
def input_binary(filename):
    clibrebound.input_binary(c_char_p(filename))
    

# Integration
def step():
    clibrebound.rebound_step()

def integrate(tmax,exactFinishTime=1,keepSynchronized=0):
    global tmpdir
    if _integrator_package == "SWIFTER":
        _particles = get_particles()
        oldwd = os.getcwd()
        paramin = """!
! Parameter file for the CHO run of the 4 giant planets and Pluto.
!
!NPLMAX         -1                 ! not used
!NTPMAX         -1                 ! not used
T0             %.16e
TSTOP          %.16e               ! length of simulation 
DT             %.16e               ! stepsize
PL_IN          pl.in
TP_IN          tp.in
IN_TYPE        ASCII
ISTEP_OUT      %d                   ! output every 10K years
BIN_OUT        bin.dat
OUT_TYPE       REAL8                ! 4-byte XDR formatted output
OUT_FORM       XV                  ! cartesian output 
OUT_STAT       NEW
ISTEP_DUMP     1000000000        ! system dump also every 10K years
J2             0.0E0               ! no J2 term
J4             0.0E0               ! no J4 term
CHK_CLOSE      no                  ! don't check for planetary close encounters
CHK_RMIN       -1.0                ! don't check for close solar encounters
CHK_RMAX       1000.0              ! discard outside of 1000 AU
CHK_EJECT      -1.0                ! ignore this check
CHK_QMIN       -1.0                ! ignore this check
!CHK_QMIN_COORD HELIO               ! commented out here
!CHK_QMIN_RANGE 1.0 1000.0          ! commented out here
ENC_OUT        enc.dat
EXTRA_FORCE    no                  ! no extra user-defined forces
BIG_DISCARD    yes                 ! output all planets if anything discarded
RHILL_PRESENT  no                  ! no Hill's sphere radii in input file
""" % ( get_t(), tmax, get_dt(),max(1,int((tmax-get_t())/get_dt())))

        if not tmpdir:
        # first call
            tmpdir = tempfile.mkdtemp()
            for f in ["swifter_whm", "swifter_tu4","swifter_symba","swifter_helio"]:
                os.symlink(oldwd+"/../../others/swifter/bin/"+f,tmpdir+"/"+f)
            os.symlink(oldwd+"/../../others/swifter/bin/tool_follow",tmpdir+"/tool_follow")
        os.chdir(tmpdir)
        smallin = """0\n"""
        with open("tp.in", "w") as f:
            f.write(smallin)
        bigin = """ %d
""" %(get_N())
        _G = get_G() 
        for i in range(0,get_N()):
            bigin += """%d    %.18e 
%.18e %.18e %.18e
%.18e %.18e %.18e
""" %(i+1, _G*_particles[i].m,_particles[i].x, _particles[i].y, _particles[i].z, _particles[i].vx, _particles[i].vy, _particles[i].vz)
        with open("pl.in", "w") as f:
            f.write(bigin)
        with open("param.in", "w") as f:
            f.write(paramin)
        #with open("param.dmp", "w") as f:
        #    f.write(paramin)
        starttime = time.time()    
        if _integrator_debug.lower() == "swifter-whm":
            os.system("echo param.in | ./swifter_whm > /dev/null")
        elif _integrator_debug.lower() == "swifter-symba":
            os.system("echo param.in | ./swifter_symba > /dev/null")
        elif _integrator_debug.lower() == "swifter-helio":
            os.system("echo param.in | ./swifter_helio > /dev/null")
        elif _integrator_debug.lower() == "swifter-tu4":
            os.system("echo param.in | ./swifter_tu4 > /dev/null")
        else:
            print("Integrator not found! %s"%_integrator_debug)
        endtime = time.time()    
        c_double.in_dll(clibrebound,"timing").value = endtime-starttime
    
        for i in range(1,get_N()):
            with open("outputparams.txt", "w") as f:
                f.write("dump_param1.dat\n")
                f.write("{0}\n".format(i+1))
                f.write("1\n")
            try:
                os.system("./tool_follow < outputparams.txt > /dev/null")
                with open("follow.out", "r") as f:
                    lines = f.readlines()
                    line = lines[-1].split()
                    t = float(line[0].strip())
                    set_t(t)
                    _particles[i].x = float(line[2].strip())
                    _particles[i].y = float(line[3].strip())
                    _particles[i].z = float(line[4].strip())
                    _particles[i].vx = float(line[5].strip())
                    _particles[i].vy = float(line[6].strip())
                    _particles[i].vz = float(line[7].strip())
            except:
                print("Something went wrong. Ignoring it for now. (%s)"%_integrator_debug)
                pass
        os.system("rm bin.dat")
        os.chdir(oldwd)
    if _integrator_package == "MERCURY":
        k = 0.01720209895    
        facTime = math.sqrt(get_G())/k
        _particles = get_particles()
        oldwd = os.getcwd()
        paramin = """)O+_06 Integration parameters  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
) Important integration parameters:
)---------------------------------------------------------------------
 algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = mvs
 start time (days)= 0.
 stop time (days) = %.16e
 output interval (days) = 1000000.1d0
 timestep (days) = %.16e
 accuracy parameter=1.d-12
)---------------------------------------------------------------------
) Integration options:
)---------------------------------------------------------------------
 stop integration after a close encounter = no
 allow collisions to occur = no
 include collisional fragmentation = no
 express time in days or years = days
 express time relative to integration start time = yes
 output precision = high
 < not used at present >
 include relativity in integration= no
 include user-defined force = no
)---------------------------------------------------------------------
) These parameters do not need to be adjusted often:
)---------------------------------------------------------------------
 ejection distance (AU)= 100
 radius of central body (AU) = 0.005
 central mass (solar) = %.16e
 central J2 = 0
 central J4 = 0
 central J6 = 0
 < not used at present >
 < not used at present >
 Hybrid integrator changeover (Hill radii) = 3.
 number of timesteps between data dumps = 50000000
 number of timesteps between periodic effects = 100000000
""" % ( tmax*facTime, get_dt()*facTime,_particles[0].m)
        if not tmpdir:
            # first call
            tmpdir = tempfile.mkdtemp()
            for f in ["mercury", "message.in","files.in"]:
                os.symlink(oldwd+"/../../others/mercury6/"+f,tmpdir+"/"+f)
            os.chdir(tmpdir)
            smallin = """)O+_06 Small-body initial data  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
 style (Cartesian, Asteroidal, Cometary) = Ast
)---------------------------------------------------------------------
"""
            with open("small.in", "w") as f:
                f.write(smallin)
            bigin = """)O+_06 Big-body initial data  (WARNING: Do not delete this line!!)
) Lines beginning with `)' are ignored.
)---------------------------------------------------------------------
 style (Cartesian, Asteroidal, Cometary) = Cartesian
 epoch (in days) = %.16e
)---------------------------------------------------------------------
""" %(get_t()*facTime)
            for i in range(1,get_N()):
                bigin += """PART%03d    m=%.18e r=20.D0 d=5.43
 %.18e %.18e %.18e
 %.18e %.18e %.18e
  0. 0. 0.
""" %(i, _particles[i].m, _particles[i].x, _particles[i].y, _particles[i].z, _particles[i].vx/facTime, _particles[i].vy/facTime, _particles[i].vz/facTime)
            with open("big.in", "w") as f:
                f.write(bigin)
            with open("param.in", "w") as f:
                f.write(paramin)
        else:
            # Not first call
            os.chdir(tmpdir)
            with open("param.dmp", "w") as f:
                f.write(paramin)
        starttime = time.time()    
        #os.system("./mercury ")
        os.system("./mercury >/dev/null")
        endtime = time.time()    
        c_double.in_dll(clibrebound,"timing").value = endtime-starttime
        with open("big.dmp", "r") as f:
            lines = f.readlines()
            t= float(lines[4].split("=")[1].strip())
            set_t(t/facTime)
            j = 1
            for i in range(6,len(lines),4):
                pos = lines[i+1].split()
                _particles[j].x = float(pos[0])
                _particles[j].y = float(pos[1])
                _particles[j].z = float(pos[2])
                vel = lines[i+2].split()
                _particles[j].vx = float(vel[0])*facTime
                _particles[j].vy = float(vel[1])*facTime
                _particles[j].vz = float(vel[2])*facTime
                j += 1

        os.chdir(oldwd)
    if _integrator_package =="REBOUND":
        clibrebound.integrate(c_double(tmax),c_int(exactFinishTime),c_int(keepSynchronized))

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


def get_com(last=None):
    """Returns the center of momentum for all particles in the simulation"""
    m = 0.
    x = 0.
    y = 0.
    z = 0.
    vx = 0.
    vy = 0.
    vz = 0.
    ps = get_particles()    # particle pointer
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

# Alias
get_center_of_momentum = get_com

# Import at the end to avoid circular dependence
from .particle import *
from . import horizons
