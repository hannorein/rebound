from ctypes import *

#Find the migration C library
pymodulespath = os.path.dirname(__file__)
try:
    libmigration = CDLL(pymodulespath + '/../shared/libmigration/libmigration.so', RTLD_GLOBAL)
except:
    print("Cannot find library 'libmigration.so'. Try typing make in rebound/shared/ and/or check path set in 'rebound/python_modules/migration.py'.")
    raise

try:
    import builtins      # if this succeeds it's python 3.x
    builtins.xrange = range
    builtins.basestring = (str,bytes)
except ImportError:
    pass                 # python 2.x

def add_migration(tau_a):
    arr = (c_double * get_N())(*tau_a)
    librebound.add_migration(byref(arr))
    librebound.set_additional_forces(CFUNCTYPE(libmigration.disk_forces))

def add_e_damping(tau_e, p=0.):
    c_double.in_dll(librebound, "e_damping_p").value = p
    arr = (c_double * get_N())(*tau_e)
    librebound.add_e_damping(byref(arr))
    
def set_e_damping(tau_e):
    arr = (c_double * get_N())(*tau_e)
    librebound.set_e_damping(byref(arr))
    
def add_i_damping(tau_i):
    arr = (c_double * get_N())(*tau_i)
    librebound.add_i_damping(byref(arr))
    
def add_peri_precession(gamma, Rc, podot):
    librebound.add_peri_precession(c_double(gamma), c_double(Rc), c_double(podot))

