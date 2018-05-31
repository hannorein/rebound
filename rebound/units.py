# Unit constants
import math

# All units entered in SI (kg, m, s)
G_SI = 6.67408e-11
times_SI = {'s':1.,
    'hr':3600.,
    'day': 86400.,
    'days': 86400.,
    'd': 86400.,
    'yr':31557600., # Julian year (exact)
    'yrs':31557600., 
    'jyr':31557600.,
    'sidereal_yr':31558149.7635,
    'yr2pi':math.sqrt(149597870700.**3/1.3271244004193938e20), # chosen to make G=1
    'kyr':31557600.*1.e3,
    'myr':31557600.*1.e6,
    'gyr':31557600.*1.e9}
lengths_SI =  {'m':1.,
    'cm':0.01,
    'km':1000.,
    'au':149597870700.,
    'aus':149597870700.
    }

    #What we measure accurately is GM, so set mass units such that G*M gives the value of GM in horizons.py (in the list at the end of horizons.py, the NAIF codes ending in 99 refer to the planets, single digits to the total mass of the planet plus its moons).  Have to multiply by 10**9 since that list has G in kg^-1km^3/s^2 and we use SI.

masses_SI = {'kg':1.,
    'g':1.0e-3,
    'gram':1.0e-3,
    'msun':1.3271244004193938E+11/G_SI*10**9,
    'msolar':1.3271244004193938E+11/G_SI*10**9,
    'mmercury':2.2031780000000021E+04/G_SI*10**9,
    'mvenus':3.2485859200000006E+05/G_SI*10**9,
    'mearth':3.9860043543609598E+05/G_SI*10**9,
    'mmars':4.282837362069909E+04/G_SI*10**9,
    'mjupiter':1.266865349218008E+08/G_SI*10**9,
    'msaturn':3.793120749865224E+07/G_SI*10**9,
    'muranus':5.793951322279009E+06/G_SI*10**9,
    'mneptune':6.835099502439672E+06/G_SI*10**9,
    'mpluto':8.696138177608748E+02/G_SI*10**9}


def units_convert_particle(p, old_l, old_t, old_m, new_l, new_t, new_m):
    p.m = convert_mass(p.m, old_m, new_m)
    p.x = convert_length(p.x, old_l, new_l) 
    p.y = convert_length(p.y, old_l, new_l)
    p.z = convert_length(p.z, old_l, new_l)
    p.vx = convert_vel(p.vx, old_l, old_t, new_l, new_t)
    p.vy = convert_vel(p.vy, old_l, old_t, new_l, new_t)
    p.vz = convert_vel(p.vz, old_l, old_t, new_l, new_t)
    p.ax = convert_acc(p.ax, old_l, old_t, new_l, new_t)
    p.ay = convert_acc(p.ay, old_l, old_t, new_l, new_t)
    p.az = convert_acc(p.az, old_l, old_t, new_l, new_t)
    return p

def convert_mass(mass, old_m, new_m):
    return mass*masses_SI[old_m]/masses_SI[new_m]

def convert_length(length, old_l, new_l):
    return length*lengths_SI[old_l]/lengths_SI[new_l]

def convert_vel(vel, old_l, old_t, new_l, new_t):
    in_SI=vel*lengths_SI[old_l]/times_SI[old_t]
    return in_SI*times_SI[new_t]/lengths_SI[new_l]

def convert_acc(acc, old_l, old_t, new_l, new_t):
    in_SI=acc*lengths_SI[old_l]/times_SI[old_t]**2
    return in_SI*times_SI[new_t]**2/lengths_SI[new_l]

def convert_G(new_l, new_t, new_m):
    return G_SI*masses_SI[new_m]*times_SI[new_t]**2/lengths_SI[new_l]**3
       
def check_units(newunits):   
    if len(newunits) is not 3:
        raise Exception("Error: Need to pass exactly 3 units for length, time, and mass (any order), see ipython_examples/Units.ipynb")
    
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
        raise Exception("Error: Need to assign rebound.units a tuple consisting of 3 units for length, time, and mass (any order).  See ipython/examples/Units.ipynb.  If you passed such a tuple, at least one of your units isn't in our list.  Please update the dictionaries at the top of rebound/rebound/units.py and send a pull request!")

    return (l_unit, t_unit, m_unit)
