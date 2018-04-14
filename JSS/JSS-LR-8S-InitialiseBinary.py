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
#				solar system planets. 
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
					   "Europa", "Ganymede", "Callisto", "Sun", "Mercury", 
					   "Venus", "Earth", "Mars", "Saturn", "Uranus", "Neptune"] 
		  
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
	model.particles[i].hash=labels[i]

# Correct mass of bodies that do not have mass information in Horizons.
solar_mass = 1.988435e+30 					# Solar Mass in kg

model.particles['Metis'].m = 3.6e16/solar_mass
model.particles['Adrastea'].m = 2e15/solar_mass
model.particles['Amalthea'].m = 7.77e17/solar_mass
model.particles['Thebe'].m = 7.17e18/solar_mass

# Calculate orbits of the bodies. Convert Cartesian coordinates to Keplerian 
# that coordinats and moves model to barycentre (the center of momentum frame
# where the center of mass is at the origin and does not move).
model.move_to_com()
model.calculate_orbits()

# Save initial conditions to binary file
model.save(BinFile)

# output all particles and integrator settings to logfile
model.status()
print('Units:\t\t{}'.format(model.units)) 
print("G:\t\t\{}.".format(model.G))
print("sim t:\t\t{} = Julian day {}".format(model.t, model.t+TimeOfEphemeris))
print("v_earth:\t{} AU/yr (Earth's orbit equals to 2pi AU)\n-------------------------------------".format(math.sqrt(model.particles["Earth"].vx**2 + model.particles["Earth"].vy**2 + model.particles["Earth"].vz**2)))


