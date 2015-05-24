# -*- coding: utf-8 -*-
"""
Functions for converting units
"""

units = {'m':1.,
        'km':1000.,
        'AU':1.5e11,
        \
        's':1.,
        'hr':3600.,
        'yr':3.15e7,
        \
        'kg':1.,
        'Msun':1.99e30
        }

def convert(val, oldunits, newunits):
    return val*units[oldunits]/units[newunits]

def convert_vel(val, oldlength, oldtime, newlength, newtime):
    in_SI_units=val*units[oldlength]/units[oldtime]
    return in_SI_units*units[newtime]/units[newlength]

def convert_G(val, oldlength, oldtime, oldmass, newlength, newtime, newmass):
    in_SI_units=val*units[oldlength]**3/units[oldmass]/units[oldtime]**2
    return in_SI_units*units[newmass]*units[newtime]**2/units[newlength]**3

