# -*- coding: utf-8 -*-
"""
Functions for converting units
"""

import rebound

units = [None, None, None]

unit_dic = {'m':1.,
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
    return val*unit_dic[oldunits]/unit_dic[newunits]

def convert_vel(val, oldlength, oldtime, newlength, newtime):
    in_SI_units=val*unit_dic[oldlength]/unit_dic[oldtime]
    return in_SI_units*unit_dic[newtime]/unit_dic[newlength]

def convert_G(val, oldlength, oldtime, oldmass, newlength, newtime, newmass):
    in_SI_units=val*unit_dic[oldlength]**3/unit_dic[oldmass]/unit_dic[oldtime]**2
    return in_SI_units*unit_dic[newmass]*unit_dic[newtime]**2/unit_dic[newlength]**3

def set_sim_units(length, time, mass)
def _units(length, time, mass):
    if sim_unit_dic[0] is None:
        if rebound.particles is not []:
            print("Error:  Have to set unit_dic BEFORE calling set_sim_unit_dic.  See unit_dic.ipynb in python_tutorials.")
            exit()

        sim_unit_dic = length, time, mass
        print(sim_unit_dic)
            
