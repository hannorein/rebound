from ctypes import *
import math
import os
import tempfile
import shutil
import time
import ctypes.util
import pkg_resources

pymodulespath = os.path.dirname(__file__)
#Find the rebound C library
try:
    librebound = CDLL(pymodulespath+"/../librebound.so", RTLD_GLOBAL)
except:
    try:
        librebound = CDLL(pymodulespath + '/../shared/librebound/librebound.so', RTLD_GLOBAL)
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
    return c_char_p.in_dll(librebound, "build_str").value

class Orbit():
    """Defines the same data structure as in tools.h"""
    def __init__(self):
        self.a      =   None    # semimajor axis
        self.r      =   None    # radial distance from reference
        self.h      =   None    # angular momentum
        self.P      =   None    # orbital period
        self.l      =   None    # mean longitude = Omega + omega + M
        self.e      =   None    # eccentricity
        self.inc    =   None    # inclination
        self.Omega  =   None    # longitude of ascending node
        self.omega  =   None    # argument of perihelion
        self.f      =   None    # true anomaly

    def __str__(self):
        return "<rebound.Orbit instance, a=%s e=%s>"%(str(self.a),str(self.e))


class Particle(Structure):
    """A particle datastructure. Same as defined in particle.h"""
    _fields_ = [("x", c_double),
                ("y", c_double),
                ("z", c_double),
                ("vx", c_double),
                ("vy", c_double),
                ("vz", c_double),
                ("ax", c_double),
                ("ay", c_double),
                ("az", c_double),
                ("m", c_double) ]
    def __str__(self):
        return "<rebound.Particle object, m=%s x=%s y=%s z=%s vx=%s vy=%s vz=%s>"%(self.m,self.x,self.y,self.z,self.vx,self.vy,self.vz)
    
    def __init__(self, particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None):   
        if particle is not None:
            raise ValueError("Cannot initialise particle from other particles.")
        cart = [x,y,z,vx,vy,vz]
        orbi = [primary,a,anom,e,omega,inc,Omega,MEAN]
        if m is None:   #default value for mass
            m = 0.
        if notNone(cart) and notNone(orbi):
                raise ValueError("You cannot pass cartesian coordinates and orbital elements at the same time.")
        if notNone(orbi):
            if primary is None:
                primary = get_center_of_momentum()
            if a is None:
                raise ValueError("You need to pass a semi major axis to initialize the particle using orbital elements.")
            if anom is None:
                anom = 0.
            if e is None:
                e = 0.
            if omega is None:
                omega = 0.
            if inc is None:
                inc = 0.
            if Omega is None:
                Omega = 0.
            if MEAN is None:
                MEAN = False
            self.set_orbit(m=m,primary=primary,a=a,anom=anom,e=e,omega=omega,inc=inc,Omega=Omega,MEAN=MEAN)
        else:
            if x is None:
                x = 0.
            if y is None:
                y = 0.
            if z is None:
                z = 0.
            if vx is None:
                vx = 0.
            if vy is None:
                vy = 0.
            if vz is None:
                vz = 0.
            self.m = m
            self.x = x
            self.y = y
            self.z = z
            self.vx = vx
            self.vy = vy
            self.vz = vz

    def set_orbit(self,
                    m,          # mass
                    primary,    # central body (rebound.Particle object)
                    a,          # semimajor axis
                    anom=0.,    # anomaly
                    e=0.,       # eccentricity
                    omega=0.,   # argument of pericenter
                    inc=0.,     # inclination
                    Omega=0.,   # longitude of ascending node
                    MEAN=False):    # mean anomaly
        """ Initialises a particle structure with the passed set of
            orbital elements. Mass (m), primary and 'a' are required (see Parameters
            below, and any orbital mechanics text, e.g., Murray & Dermott
            Solar System Dynamics for definitions). Other values default to zero.
            All angles should be passed in radians. Units are set by the
            gravitational constant G (default = 1.). If MEAN is set to True, anom is
            taken as the mean anomaly, rather than the true anomaly.
            
            Usage
            _____
            TODO: UPDATE
            primary = rebound.Particle(m=1.) # particle with unit mass at origin & v=0
            
            # test particle (m=0) with specified elements using mean anomaly
            p = kepler_particle(m=0.,primary=primary,a=2.5, anom=math.pi/2,e=0.3,
            omega=math.pi/6,inc=math.pi/3,Omega=0.,MEAN=True)
            
            # m=0.1 particle on circular orbit math.pi/4 from x axis in xy plane
            p = kepler_particle(0.1,primary,2.5,math.pi/4)
            
            Parameters
            __________
            m       : (float)            Mass of the particle
            primary : (rebound.Particle) Particle structure for the central body
            a       : (float)            Semimajor axis
            anom    : (float)            True anomaly (default).
            Mean anomaly if MEAN is set to True
            e       : (float)            Eccentricity
            omega   : (float)            Argument of pericenter
            inc     : (float)            Inclination (to xy plane)
            Omega   : (float)            Longitude of the ascending node
            MEAN    : (boolean)          If False (default), anom = true anomaly
            If True, anom = mean anomaly
            
            Returns
            _______
            A rebound.Particle structure initialized with the given orbital parameters
            """
        
        self.m = m

        if not(0.<=inc<=math.pi): raise ValueError('inc must be in range [0,pi]')
        if e>1.:
            if math.fabs(anom)>math.acos(-1./e): raise ValueError('hyperbolic orbit with anomaly larger than angle of asymptotes')
        
        if MEAN is True: # need to calculate f
            E = eccentricAnomaly(e,anom)
            if e>1.:
                print("not working yet")
                exit(1)
                f = 2.*math.atan(math.sqrt(-(1.+ e)/(1. - e))*math.tanh(0.5*E))
            else:
                f = mod2pi(2.*math.atan(math.sqrt((1.+ e)/(1. - e))*math.tan(0.5*E)))
        else:
            f = anom
        
        cO = math.cos(Omega)
        sO = math.sin(Omega)
        co = math.cos(omega)
        so = math.sin(omega)
        cf = math.cos(f)
        sf = math.sin(f)
        ci = math.cos(inc)
        si = math.sin(inc)
        
        r = a*(1.-e**2)/(1.+e*cf)
        
        # Murray & Dermott Eq. 2.122
        self.x  = primary.x + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci)
        self.y  = primary.y + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci)
        self.z  = primary.z + r*(so*cf+co*sf)*si
        
        n = math.sqrt(get_G()*(primary.m+m)/(a**3))
        if e>1.:
            v0 = n*a/math.sqrt(-(1.-e**2))
        else:
            v0 = n*a/math.sqrt(1.-e**2)
        
        # Murray & Dermott Eq. 2.36 after applying the 3 rotation matrices from Sec. 2.8 to the velocities in the orbital plane
        self.vx = primary.vx + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO))
        self.vy = primary.vy + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO))
        self.vz = primary.vz + v0*((e+cf)*co*si - sf*si*so)


    def get_orbit(self, primary=None, verbose=False):
        """ Returns a rebound.Orbit object with the keplerian orbital elements
            corresponding to the particle around the central body primary
            (rebound.Particle). Edge cases will return values set to None. If
            verbose is set to True (default=False), error messages are printed
            when a breakout condition is met.
            
            Usage
            _____
            TODO: Update!
            orbit = p2orbit(p,primary)
            print(orbit.e) # gives the eccentricity
            
            orbit = p2orbit(p,primary,verbose=True) # will print out error msgs
            
            Parameters
            __________
            self     : (rebound.Particle) particle for which orbital elements are sought
            primary  : (rebound.Particle) central body
            verbose  : (boolean)          If set to True, will print out error msgs
            
            Returns
            _______
            A rebound.Orbit object (with member variables for the orbital elements)
            """
        if primary is None:
            primary = get_particles()[0]

        o = Orbit()
        if primary.m <= TINY:
            if verbose is True:
                print("Star has no mass.")
            return o                            # all values set to None
        
        dx = self.x - primary.x
        dy = self.y - primary.y
        dz = self.z - primary.z
        o.r = math.sqrt ( dx*dx + dy*dy + dz*dz )
        if o.r <= TINY:
            if verbose is True:
                print('Particle and primary positions are the same.')
            return o
        
        dvx = self.vx - primary.vx
        dvy = self.vy - primary.vy
        dvz = self.vz - primary.vz
        v = math.sqrt ( dvx*dvx + dvy*dvy + dvz*dvz )
        
        mu = get_G()*(self.m+primary.m)
        o.a = -mu/( v*v - 2.*mu/o.r )               # semi major axis
        
        h0 = (dy*dvz - dz*dvy)                      # angular momentum vector
        h1 = (dz*dvx - dx*dvz)
        h2 = (dx*dvy - dy*dvx)
        o.h = math.sqrt ( h0*h0 + h1*h1 + h2*h2 )   # abs value of angular momentum
        if o.h/(o.r*v) <= MIN_REL_ERROR:
            if verbose is True:
                print('Particle orbit is radial.')
            return o
        
        vr = (dx*dvx + dy*dvy + dz*dvz)/o.r
        e0 = 1./mu*( (v*v-mu/o.r)*dx - o.r*vr*dvx )
        e1 = 1./mu*( (v*v-mu/o.r)*dy - o.r*vr*dvy )
        e2 = 1./mu*( (v*v-mu/o.r)*dz - o.r*vr*dvz )
        o.e = math.sqrt( e0*e0 + e1*e1 + e2*e2 )   # eccentricity
        
        o.P = math.copysign(2.*math.pi*math.sqrt( math.fabs(o.a*o.a*o.a/mu) ), o.a)  # period
        o.inc = math.acos( h2/o.h )               # inclination (wrt xy-plane)
        # if pi/2<i<pi it's retrograde
        n0 = -h1                                    # node vector
        n1 =  h0                                    # in xy plane => no z component
        n = math.sqrt( n0*n0 + n1*n1 )
        er = dx*e0 + dy*e1 + dz*e2
        if n/(o.r*v)<=MIN_REL_ERROR or o.inc<=MIN_REL_ERROR:# we are in the xy plane
            o.Omega=0.
            if o.e <= MIN_REL_ERROR:              # omega not defined for circular orbit
                o.omega = 0.
            else:
                if e1>=0.:
                    o.omega=math.acos(e0/o.e)
                else:
                    o.omega = 2.*math.pi-math.acos(e0/o.e)
        else:
            if o.e <= MIN_REL_ERROR:
                o.omega = 0.
            else:
                if e2>=0.:                        # omega=0 if perictr at asc node
                    o.omega=math.acos(( n0*e0 + n1*e1 )/(n*o.e))
                else:
                    o.omega=2.*math.pi-math.acos(( n0*e0 + n1*e1 )/(n*o.e))
            if n1>=0.:
                o.Omega = math.acos(n0/n)
            else:
                o.Omega=2.*math.pi-math.acos(n0/n)# Omega=longitude of asc node
        # taken in xy plane from x axis
        
        if o.e<=MIN_REL_ERROR:                           # circular orbit
            o.f=0.                                  # f has no meaning
            o.l=0.
        else:
            cosf = er/(o.e*o.r)
            cosea = (1.-o.r/o.a)/o.e
            
            if -1.<=cosf and cosf<=1.:                       # failsafe
                o.f = math.acos(cosf)
            else:
                o.f = math.pi/2.*(1.-cosf)
            
            if -1.<=cosea and cosea<=1.:
                ea  = math.acos(cosea)
            else:
                ea = math.pi/2.*(1.-cosea)
            
            if vr<0.:
                o.f=2.*math.pi-o.f
                ea =2.*math.pi-ea
            
            o.l = ea -o.e*math.sin(ea) + o.omega+ o.Omega  # mean longitude
        
        return o

def status():
    """ Returns a string with a summary of the current status 
        of the simulation
        """
    s= ""
    N = get_N()
    s += "---------------------------------\n"
    s += "Rebound version: "+ pkg_resources.require("rebound")[0].version +"\n"
    s += get_build_str() + "\n"
    s += "Number of particles: \t%d\n" %N       
    s += "Simulation time:     \t%f\n" %get_t()
    if N>0:
        s += "---------------------------------\n"
        p = get_particles()
        for i in range(N):
            s += str(p[i]) + "\n"
    s += "---------------------------------"
    return s

# Set function pointer for additional forces

AFF = CFUNCTYPE(None)
fp = None
def set_additional_forces(func):
    global fp  # keep references
    fp = AFF(func)
    librebound.set_additional_forces(fp)

# Setter/getter of parameters and constants
def set_G(G):
    c_double.in_dll(librebound, "G").value = G

def get_G():
    return c_double.in_dll(librebound, "G").value

def set_dt(dt):
    c_double.in_dll(librebound, "dt").value = dt

def get_dt():
    return c_double.in_dll(librebound, "dt").value

def set_t(t):
    c_double.in_dll(librebound, "t").value = t

def set_min_dt(t):
    c_double.in_dll(librebound, "integrator_ias15_min_dt").value = t

def get_t():
    return c_double.in_dll(librebound, "t").value

def init_megno(delta):
    librebound.tools_megno_init(c_double(delta))

def get_megno():
    librebound.tools_megno.restype = c_double
    return librebound.tools_megno()

def get_lyapunov():
    librebound.tools_lyapunov.restype = c_double
    return librebound.tools_lyapunov()

def get_N():
    return c_int.in_dll(librebound,"N").value 

def get_N_megno():
    return c_int.in_dll(librebound,"N_megno").value 

def get_iter():
    return c_int.in_dll(librebound,"iter").value 

def get_timing():
    return c_double.in_dll(librebound,"timing").value 

# Setter functions, particle data
def set_particles(_particles):
    c_int.in_dll(librebound,"N").value = len(_particles)
    arr = (Particle * len(_particles))(*_particles)
    librebound.setp(byref(arr))

def add_particle(particle=None, m=None, x=None, y=None, z=None, vx=None, vy=None, vz=None, primary=None, a=None, anom=None, e=None, omega=None, inc=None, Omega=None, MEAN=None):   
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
            librebound.particles_add(particle)
    else:
       librebound.particles_add(particle)

# Aliases
add_particles = add_particle
add = add_particle

# Particle getter functions
def get_particle(i):
    N = get_N() 
    if i>=N:
        return None
    getp = librebound.particle_get
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
    N = c_int.in_dll(librebound,"N").value 
    getp = librebound.particles_get
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
            com = get_center_of_momentum(i)
        orbits.append(_particles[i].get_orbit(primary=com))
    return orbits

# Tools
def move_to_center_of_momentum():
    librebound.tools_move_to_center_of_momentum()

tmpdir = None
def reset():
    global tmpdir
    if tmpdir:
        shutil.rmtree(tmpdir)
        tmpdir = None
    librebound.reset()

def set_integrator_mikkola_corrector(on=0):
    c_int.in_dll(librebound, "integrator_mikkola_corrector").value = on

#REMOVE Following variables
integrator_package = "REBOUND"
integrator_debug = ""    

def set_integrator(integrator="IAS15"):
    global integrator_package
    global integrator_debug
    if isinstance(integrator, int):
        librebound.integrator_set(c_int(integrator))
        return
    if isinstance(integrator, basestring):
        integrator_debug = integrator
        integrator_package = "REBOUND"
        if integrator.lower() == "ias15":
            set_integrator(0)
            return
        if integrator[0:7].lower() == "mikkola":
            set_integrator(1)
            set_integrator_mikkola_corrector(0)
            if integrator[8:] == "cor3":
                set_integrator_mikkola_corrector(3)
            if integrator[8:] == "cor5":
                set_integrator_mikkola_corrector(5)
            if integrator[8:] == "cor7":
                set_integrator_mikkola_corrector(7)
            if integrator[8:] == "cor11":
                set_integrator_mikkola_corrector(11)
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
            integrator_package = "MERCURY"
            return
        if integrator.lower() == "swifter-whm":
            integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-symba":
            integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-helio":
            integrator_package = "SWIFTER"
            return
        if integrator.lower() == "swifter-tu4":
            integrator_package = "SWIFTER"
            return
    raise ValueError("Warning. Intergrator not found.")

def set_integrator_mikkola_persistent_particles(is_per=0):
    if isinstance(is_per, int):
        c_int.in_dll(librebound, "integrator_mikkola_persistent_particles").value = is_per
        return
    raise ValueError("Expecting integer.")

def set_integrator_mikkola_synchronize_manually(synchronize_manually=0):
    if isinstance(synchronize_manually, int):
        c_int.in_dll(librebound, "integrator_mikkola_synchronize_manually").value = synchronize_manually
        return
    raise ValueError("Expecting integer.")

def set_force_is_velocitydependent(force_is_velocitydependent=1):
    if isinstance(force_is_velocitydependent, int):
        c_int.in_dll(librebound, "integrator_force_is_velocitydependent").value = force_is_velocitydependent
        return
    raise ValueError("Expecting integer.")

# Input/Output routines

def output_binary(filename):
    librebound.output_binary(c_char_p(filename))
    
def input_binary(filename):
    librebound.input_binary(c_char_p(filename))
    

# Integration
def step():
    librebound.rebound_step()

def integrate(tmax,exactFinishTime=1,keepSynchronized=0):
    global tmpdir
    if integrator_package == "SWIFTER":
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
        if integrator_debug.lower() == "swifter-whm":
            os.system("echo param.in | ./swifter_whm > /dev/null")
        elif integrator_debug.lower() == "swifter-symba":
            os.system("echo param.in | ./swifter_symba > /dev/null")
        elif integrator_debug.lower() == "swifter-helio":
            os.system("echo param.in | ./swifter_helio > /dev/null")
        elif integrator_debug.lower() == "swifter-tu4":
            os.system("echo param.in | ./swifter_tu4 > /dev/null")
        else:
            print("Integrator not found! %s"%integrator_debug)
        endtime = time.time()    
        c_double.in_dll(librebound,"timing").value = endtime-starttime
    
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
                print("Something went wrong. Ignoring it for now. (%s)"%integrator_debug)
                pass
        os.system("rm bin.dat")
        os.chdir(oldwd)
    if integrator_package == "MERCURY":
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
        c_double.in_dll(librebound,"timing").value = endtime-starttime
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
    if integrator_package =="REBOUND":
        librebound.integrate(c_double(tmax),c_int(exactFinishTime),c_int(keepSynchronized))

TWOPI = 2.*math.pi
def mod2pi(f):
    """Returns the angle f modulo 2 pi."""
    while f<0.:
        f += TWOPI
    while f>TWOPI:
        f -= TWOPI
    return f

def notNone(a):
    """Returns True if array a contains at least one element that is not None. Returns False otherwise."""
    return a.count(None) != len(a)

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


def get_center_of_momentum(last=None):
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


