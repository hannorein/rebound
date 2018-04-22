#!/usr/bin/env python
# 
###############################################################################
# Program:		JSS-LR-8S-InCond.py
#				(Jupiter Satellite System Long Run 8 Satellites)
#
# Author:		Christopher C.E. Tylor
# 				University of Southern Queensland
# EMail:		christopher.tylor@usq.edu.au
#
# Version:		1.0
# 				14. April 2018
# Based on:		rebound by Hanno Rein http://github.com/hannorein/rebound
# Usage: 		python JSS-LR-8S-InCond.py 
#------------------------------------------------------------------------------
# Description:	This program initialises the binary file used for the analysis 
#				of the orbital evolution of Jupiter's inner 8 satellites. The 
#				file also contains the initial conditions of the Sun and all 
#				outer solar system planets. The terrestrial planets are added 
#				to the Sun's mass. 
#				The program needs only to be run once to create the file
#				JSS-LR-8S-InCond-JD2458208.bin 
#				The binary file is used as base for all simulations. The 
#				ephemeris is fixed at 30-03-2018 12:00 UTC.
###############################################################################

import rebound as rb
import datetime
import math
import re
import utils

###############################################################################
# define variables pretending to be constants
###############################################################################
MASS_SOLAR 			= addons.getMass('Sun')				# in kg
MASS_MERCURY 		= addons.getMass('Mercury')			# in kg
MASS_VENUS 			= addons.getMass('Venus')			# in kg
MASS_EARTH 			= addons.getMass('Earth barycenter')	# in kg (Earth+Moon)
MASS_MARS 			= addons.getMass('Mars barycenter')	# in kg (Mars+Moons)

###############################################################################
# Set bodies for model
###############################################################################
bodies  =  ["Jupiter", "Metis", "Adrastea", "Amalthea", "Thebe", "Io", "Europa", 
		    "Ganymede", "Callisto", "Sun", "Saturn barycenter", "Uranus barycenter", 
		    "Neptune barycenter"]

for i in range(len(bodies)):bodies[i] = str(addons.getNAIF(bodies[i]))

###############################################################################
# Set rebound integrator conditions
###############################################################################
model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")	# chosen yr instead of yr/2pi
model.integrator 	= "whfast"				# whfast integrator (no collision)
#model.integrator 	= "IAS15"				# IAS15 integrator (collision)
time_of_ephemeris	= 2458208.				# Julian Day of 30/03/2018 12:00
model.t 			= time_of_ephemeris
bin_file 			= '../data/JSS-LR-8S-InCond-JD'+str(time_of_ephemeris)+'.bin'

model.add(bodies, date='2018-03-30 12:00') 	# Julian day 2458208.

###############################################################################
# Set particle hash to address particles by their name rather than their index
###############################################################################
for i in range(0,len(model.particles),1):
	model.particles[i].hash = utils.getNAIF(int(bodies[i]))

###############################################################################
# Throw the inner planets into the Sun "Boohahaha!"
###############################################################################
model.particles['Sun'].m 		= (MASS_SOLAR+ 		\
	                               MASS_MERCURY+	\
	                               MASS_VENUS+		\
	                               MASS_EARTH+		\
	                               MASS_MARS)/		\
								   MASS_SOLAR

###############################################################################
# Save initial conditions to binary file
###############################################################################
model.status()
model.save(bin_file)
