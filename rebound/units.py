# -*- coding: utf-8 -*-
"""
Functions for converting units
"""

class Units():
    G = 6.674e-11
    units = {'m':1.,
        'km':1000.,
        'au':149597870700,
        \
        's':1.,
        'hr':3600.,
        'yr':31557600.0, # Julian year (exact)
        'kyr':31557600.0*1.e3,
        'myr':31557600.0*1.e6,
        'gyr':31557600.0*1.e9,
        \
        #What we measure accurately is GM, so set mass units such that G*M gives the value of GM in horizons.py (in the list at the end of horizons.py, the NAIF codes ending in 99 refer to the planets, single digits to the total mass of the planet plus its moons).  Have to multiply by 10**9 since that list has G in kg^-1km^3/s^2 and we use SI.
        'kg':1.,
        'msun':1.3271244004193938E+11/G*10**9,
        'mmercury':2.2031780000000021E+04/G*10**9,
        'mvenus':3.2485859200000006E+05/G*10**9,
        'mearth':3.9860043543609598E+05/G*10**9,
        'mmars':4.282837362069909E+04/G*10**9,
        'mjupiter':1.266865349218008E+08/G*10**9,
        'msaturn':3.793120749865224E+07/G*10**9,
        'muranus':5.793951322279009E+06/G*10**9,
        'mneptune':6.835099502439672E+06/G*10**9,
        'mpluto':8.696138177608748E+02/G*10**9
        
        }

    sim_units = {'length':None, 'time':None, 'mass':None}

    def __init__(self, newlength, newtime, newmass):
        import rebound
        if rebound.particles: # some particles are loaded
            print("Error:  You must initialize the units before populating the particles array.  See Units.ipynb in python_tutorials.")
            import sys
            sys.exit()

        self.sim_units['length'] = newlength.lower()
        self.sim_units['time'] = newtime.lower()
        self.sim_units['mass'] = newmass.lower()
        rebound.G = self.convert_G()

    def switch(self, newlength, newtime, newmass):
        import rebound
       
        for p in rebound.particles:
            p.m = self.convert(p.m, self.sim_units['mass'], newmass.lower())
            p.x = self.convert(p.x, self.sim_units['length'], newlength.lower())
            p.y = self.convert(p.y, self.sim_units['length'], newlength.lower())
            p.z = self.convert(p.z, self.sim_units['length'], newlength.lower())
            p.vx = self.convert_vel(p.vx, self.sim_units['length'], self.sim_units['time'], newlength.lower(), newtime.lower())
            p.vy = self.convert_vel(p.vy, self.sim_units['length'], self.sim_units['time'], newlength.lower(), newtime.lower())
            p.vz = self.convert_vel(p.vz, self.sim_units['length'], self.sim_units['time'], newlength.lower(), newtime.lower())
            p.ax = self.convert_acc(p.ax, self.sim_units['length'], self.sim_units['time'], newlength.lower(), newtime.lower())
            p.ay = self.convert_acc(p.ay, self.sim_units['length'], self.sim_units['time'], newlength.lower(), newtime.lower())
            p.az = self.convert_acc(p.az, self.sim_units['length'], self.sim_units['time'], newlength.lower(), newtime.lower())

        self.sim_units['length'] = newlength.lower()
        self.sim_units['time'] = newtime.lower()
        self.sim_units['mass'] = newmass.lower()

        rebound.G = self.convert_G()

    def convert(self, val, oldunits, newunits):
        return val*self.units[oldunits]/self.units[newunits]

    def convert_vel(self, val, oldlength, oldtime, newlength, newtime):
        in_SI_units=val*self.units[oldlength]/self.units[oldtime]
        return in_SI_units*self.units[newtime]/self.units[newlength]

    def convert_acc(self, val, oldlength, oldtime, newlength, newtime):
        in_SI_units=val*self.units[oldlength]/self.units[oldtime]**2
        return in_SI_units*self.units[newtime]**2/self.units[newlength]
    
    def convert_G(self):
        return self.G*self.units[self.sim_units['mass']]*self.units[self.sim_units['time']]**2/self.units[self.sim_units['length']]**3

