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

# Define "semi-constants" (Name is always in CAPS)
MASS_SOLAR 			= rb.units.masses_SI['msun'] 		# one Solar mass in kg
MASS_MERCURY		= rb.units.masses_SI['mmercury']	# in kg
MASS_VENUS			= rb.units.masses_SI['mvenus']		# in kg
MASS_EARTH			= rb.units.masses_SI['mearth']		# in kg
MASS_MARS			= rb.units.masses_SI['mmars']		# in kg

# set bodies for model (505 is Amalthea) and set rebound integrator conditions
labels 				= ["Jupiter", "Metis", "Adrastea", "505", "Thebe", "Io", 
					   "Europa", "Ganymede", "Callisto", "Sun","Saturn", 
					   "Uranus", "Neptune"] 
		  
model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")	# chosen yr instead of yr/2pi
model.integrator 	= "whfast"				# whfast integrator (no collision)
#model.integrator 	= "IAS15"				# IAS15 integrator (collision)
time_of_ephemeris	= 2458208.				# Julian Day of 30/03/2018 12:00
bin_file 			= '../data/JSS-LR-8S-InCond-JD'+str(time_of_ephemeris)+'.bin'

model.add(labels, date='2018-03-30 12:00') 	# Julian day 2458208.

# Set particle hash to address particles by their name rather than their index
labels[labels.index("505")] = "Amalthea" 	# Correct label for Amalthea

for i in range(0,len(model.particles),1):
	model.particles[i].hash = labels[i]

# Correct mass of bodies that do not have mass information in Horizons.


model.particles['Metis'].m 		= 3.6e16/MASS_SOLAR
model.particles['Adrastea'].m 	= 2e15/MASS_SOLAR
model.particles['Amalthea'].m 	= 7.77e17/MASS_SOLAR
model.particles['Thebe'].m 		= 7.17e18/MASS_SOLAR

# Throw the inner planets into the Sun "Boohahaha!"
model.particles['Sun'].m 		= 1.000005976997955015050933980092
print(model.particles['Sun'].m)
model.particles['Sun'].m 		= MASS_SOLAR/(MASS_SOLAR+ 	\
	                                          MASS_MERCURY+	\
	                                          MASS_VENUS+	\
	                                          MASS_EARTH+	\
	                                          MASS_MARS)

print(model.particles['Sun'].m)
# Save initial conditions to binary file
print('---------------------------------\nmodel.G:\t\t{}'.format(model.G))
model.status()
model.save(bin_file)
