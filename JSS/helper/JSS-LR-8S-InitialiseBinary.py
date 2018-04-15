###############################################################################
#
# Program:		JSS-LR-8S-InitialiseBinary.py
#				(Jupiter Satellite System Long Run 12 Satellites)
# Author:		Christopher C.E. Tylor
# 				University of Southern Queensland
# EMail:		christopher.tylor@usq.edu.au
# Version:		1.0
# 				14. April 2018
# Based on:		rebound by Hanno Rein http://github.com/hannorein/rebound
# Usage: python JSS-LR-8S-InitialiseBinary.py > JSS-LR-8S-InitialiseBinary.log
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

# set bodies for model (505 is Amalthea) and set rebound integrator conditions
labels 				= ["Jupiter", "Metis", "Adrastea", "505", "Thebe", "Io", 
					   "Europa", "Ganymede", "Callisto", "Sun","Saturn", 
					   "Uranus", "Neptune"] 
		  
model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")	# chosen yr instead of yr/2pi
model.integrator 	= "whfast"				# whfast integrator (no collision)
#model.integrator 	= "IAS15"				# IAS15 integrator (collision)
TimeOfEphemeris		= 2458208				# Julian Day of 30/03/2018 12:00
BinFile 			= 'JSS-LR-8S-InCond-JD'+str(TimeOfEphemeris)+'.bin'

model.add(labels, date='2018-03-30 12:00') 	# Julian day 2458208

# Set particle hash to address particles by their name rather than their index
labels[labels.index("505")] = "Amalthea" 	# Correct label for Amalthea

for i in range(0,len(model.particles),1):
	model.particles[i].hash = labels[i]

# Correct mass of bodies that do not have mass information in Horizons.
solar_mass = 1.988435e+30 	# one Solar mass in kg

model.particles['Metis'].m 		= 3.6e16/solar_mass
model.particles['Adrastea'].m 	= 2e15/solar_mass
model.particles['Amalthea'].m 	= 7.77e17/solar_mass
model.particles['Thebe'].m 		= 7.17e18/solar_mass

# Throw the inner planets into the Sun "Boohahaha!"
# Mass Sun:						  1.000000000000000000000000000000
# Mass Mercury:					  0.000000166011415305434853004009
# Mass Venus:					  0.000002447838287784771539704097
# Mass Earth:					  0.000003040432648022641606594830
# Mass Mars:					  0.000000322715603755499713418687
# ----------------------------------------------------------------
model.particles['Sun'].m 		= 1.000005976997955015050933980092

# Save initial conditions to binary file
model.save(BinFile)

