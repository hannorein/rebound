import rebound as rb
import datetime

# set bodies for model (505 is Amalthea) and set rebound integrator conditions
labels 				= ["Sun", "Mercury", "Venus", "Earth", "Mars"] 
		  
model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")	# chosen yr instead of yr/2pill
model.integrator 	= "whfast"				# whfast integrator (no collision)

model.add(labels, date='2018-03-30 12:00') 	# Julian day 2458208

# Set particle hash to address particles by their name rather than their index

for i in range(0,len(model.particles),1):
	model.particles[i].hash=labels[i]

print('\nMass Sun:\t{0:.30f}\nMass Mercury:\t{1:.30f}\nMass Venus:\t{2:.30f}\nMass Earth:\t{3:.30f}\nMass Mars:\t{4:.30f}'.format(model.particles['Sun'].m,model.particles['Mercury'].m,model.particles['Venus'].m,model.particles['Earth'].m,model.particles['Mars'].m,))
print('---------------------------------------------------------\nTotal:\t\t{0:.30f}'.format(model.particles['Sun'].m+model.particles['Mercury'].m+model.particles['Venus'].m+model.particles['Earth'].m+model.particles['Mars'].m,))

