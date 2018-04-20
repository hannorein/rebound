import rebound as rb
import datetime
import math
import re
import utils

bodies  =  ["Jupiter", "Metis"]

for i in range(len(bodies)):bodies[i] = str(utils.getNAIF(bodies[i]))
model = rb.Simulation()
model.add(bodies, date='2018-03-30 12:00') 	# Julian day 2458208.
for i in range(len(model.particles)):
	model.particles[i].hash = utils.getNAIF(int(bodies[i]))
	print(model.particles[i].m)
