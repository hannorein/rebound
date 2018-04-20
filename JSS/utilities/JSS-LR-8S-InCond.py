#!/usr/bin/env python
# 
###############################################################################
#
# Program:		JSS-LR-8S-InCond.py
#				(Jupiter Satellite System Long Run 12 Satellites)
# Author:		Christopher C.E. Tylor
# 				University of Southern Queensland
# EMail:		christopher.tylor@usq.edu.au
# Version:		1.0
# 				14. April 2018
# Based on:		rebound by Hanno Rein http://github.com/hannorein/rebound
# Usage: 		python JSS-LR-8S-InCond.py > ../data/JSS-LR-8S-InCond.py.log
#------------------------------------------------------------------------------
# Description:	This program initialises the binary file used for the analysis 
#				of the orbital evolution of Jupiter's inner 8 satellites. The 
#				file also contains the initial conditions of the Sun and all 
#				outer solar system planets. The terrestrial planets are added 
#				to the Sun's mass. 
#				The program needs only to be run once to create the file
#				JSS-LR-8S-InCond-JD2458208.bin and JSS-LR-8S-InitialiseBinary.log
#				The binary file is used as base for all simulations, the log 
#				file is a human readable record of the initial conditions. 
#				The ephemeris is fixed at 30-03-2018 12:00 UTC.
#
###############################################################################

import rebound as rb
import datetime
import math
import re
import utils

# =============================================================================
# A word of warning. NATIF makes a clear distinction between "Jupiter" which 
# is code 599 and "Jupiter barycentre" which is code 5. This is also true for
# other planets with a satellite system (i.e. 3 for Earth barycenter and 399 
# for Earth). Rebound's horizons.py module pulls always the single digit NAIF 
# code (i.e. 5 for Jupiter) unless one specifies and states the three digit 
# code ending in 99 (i.e. 599 for Jupiter) in the add call.
#
# If a planet is pulled from Horizons database by it's name, rebound pulls 
# the mass of the planet as the total mass of the planet plus the mass of 
# it's satellites. Therefore for example: "model.add('Jupiter')" will add 
# the mass of the whole Jupiter System, while # "model.add('599')" adds 
# Jupiter without it's satellites. 
# 
# To avoid this use utils.getNAIF("Name_Of_Planet") to pull the right code.
# i.e utils.getNAIF("Jupiter") for the planet
#	  utils.getNAIF("Jupiter barycenter") for the whole Jupiter system
#  
# As Horizons does not provide masses rebound defines them in the list at the 
# end of rebounds horizons.py. The NAIF codes ending in 99 refer to the planets, 
# single digits to the total mass of the planet plus its moons. 
# =============================================================================

# define variables pretending to be constants

MASS_SOLAR 			= utils.getMass('Sun')				# in kg
MASS_MERCURY 		= utils.getMass('Mercury')			# in kg
MASS_VENUS 			= utils.getMass('Venus')			# in kg
MASS_EARTH 			= utils.getMass('Earth barycenter')	# in kg
MASS_MARS 			= utils.getMass('Mars barycenter')	# in kg

# Set bodies for model

bodies  =  ["Jupiter", "Metis", "Adrastea", "Amalthea", "Thebe", "Io", "Europa", 
		    "Ganymede", "Callisto", "Sun", "Saturn barycenter", "Uranus barycenter", 
		    "Neptune barycenter"]

for i in range(0,len(bodies),1):bodies[i] = str(utils.getNAIF(bodies[i]))

# Set rebound integrator conditions

model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")	# chosen yr instead of yr/2pi
model.integrator 	= "whfast"				# whfast integrator (no collision)
#model.integrator 	= "IAS15"				# IAS15 integrator (collision)
model.t 			= 2458208.	
time_of_ephemeris	= 2458208.				# Julian Day of 30/03/2018 12:00
bin_file 			= '../data/JSS-LR-8S-InCond-JD'+str(time_of_ephemeris)+'.bin'

model.add(bodies, date='2018-03-30 12:00') 	# Julian day 2458208.

# Set particle hash to address particles by their name rather than their index

for i in range(0,len(model.particles),1):
	model.particles[i].hash = utils.getNAIF(int(bodies[i]))

print(model.particles['Metis'].m)	
print(model.particles['Adrastea'].m)
print(model.particles['Amalthea'].m)
print(model.particles['Thebe'].m) 
print('-'*40)	


# Correct mass of bodies that do not have mass information in Horizons.
# =============================================================================
# 1. Metis 
#	 This moon has an irregular shape and measures 60×40×34 km across. A very 
#    rough estimate of its surface area could be placed between 5,800 and 
#    11,600 square kilometers (approx. 8,700). The bulk composition and mass 
#    of Metis are not known, but assuming that its mean density is like that 
#    of Amalthea (~0.86 g/cm^3), its mass can be estimated as ~3.6×10^16 kg. 
# 2. Adrastea
#	 has also an irregular shape and measures 20×16×14 km across.A surface area 
#    estimate would be between 840 and 1,600 (~1,200) km2. The bulk, composition 
#    and mass of Adrastea are not known, but assuming that its mean density is 
#    like that of Amalthea, around 0.86 g/cm3,its mass can be estimated at 
#    about 2×10^15 kg.
# 3. Thebe
#    is irregularly shaped, with the closest ellipsoidal approximation being 
#	 116×98×84 km. Its surface area is probably between 31,000 and 59,000 
#    (~45,000) km^2. Its bulk density and mass are not known, but assuming 
#    that its mean density is like that of Amalthea (around 0.86 g/cm^3),its 
#    mass can be estimated at roughly 4.3 × 10^17 kg.

print('Metis:\t\t{0:.15e}\t{1:.15e}'.format(model.particles['Metis'].m,3.6e16/MASS_SOLAR))	
print('Adrastea:\t{0:.15e}\t{1:.15e}'.format(model.particles['Adrastea'].m,2e15/MASS_SOLAR))
print('Amalthea\t{0:.15e}\t{1:.15e}'.format(model.particles['Amalthea'].m,2.08e18/MASS_SOLAR))
print('Thebe\t\t{0:.15e}\t{1:.15e}'.format(model.particles['Thebe'].m,7.77e17/MASS_SOLAR)) 	
print('-'*40)
print('Metis: 516\t{0:.15e}'.format(3.6e16*rb.horizons.Gkmkgs))	
print('Adrastea 515:\t{0:.15e}'.format(2e15*rb.horizons.Gkmkgs))
print('Amalthea 505\t{0:.15e}'.format(2.08e18*rb.horizons.Gkmkgs))
print('Thebe 514\t{0:.15e}'.format(7.77e17*rb.horizons.Gkmkgs)) 	

# Throw the inner planets into the Sun "Boohahaha!"
model.particles['Sun'].m 		= (MASS_SOLAR+ 		\
	                               MASS_MERCURY+	\
	                               MASS_VENUS+		\
	                               MASS_EARTH+		\
	                               MASS_MARS)/		\
								   MASS_SOLAR

# Save initial conditions to binary file
model.status()
model.save(bin_file)
