import rebound as rb
import datetime
import math


labels = ["Jupiter", "Metis", "Adrastea", "505", "Thebe", "Io", "Europa", 
		  "Ganymede", "Callisto", "Leda", "Himalia", "Lysithea", "Elara", 
		  "Sun", "Mercury", "Venus", "Earth", "Mars", "Saturn", "Uranus", 
		  "Neptune"] # 505 is Amalthea
		  

# Initialise Model 

print("Creating new simulation with the following bodies: \n{}".format(labels))

model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")
model.integrator 	= "whfast"
TimeOfEphemeris		= 2458208
BinFile 			= 'JSS-LR-12S-InCond-JD'+str(TimeOfEphemeris)+'.bin'

print('{}'.format(BinFile))

model.add(labels, date='2018-03-30 12:00') #Julian day 2458208

print('Adding hashes...')
labels[labels.index("505")] = "Amalthea" # Correct label for Amalthea
for i in range(0,len(model.particles),1):
	model.particles[i].hash=labels[i]

print('Adding mass for the bodies that have none stored ...')

solar_mass = 1.988435e+30 #kilograms

model.particles['Metis'].m = 3.6e16/solar_mass
model.particles['Adrastea'].m = 2e15/solar_mass
model.particles['Amalthea'].m = 7.77e17/solar_mass
model.particles['Thebe'].m = 7.17e18/solar_mass
model.particles['Leda'].m = 5.58e15/solar_mass
model.particles['Himalia'].m = 6.7e18/solar_mass
model.particles['Lysithea'].m = 7.77e16/solar_mass
model.particles['Elara'].m = 8.7e17/solar_mass

model.status()

print('Units:\t\t{}'.format(model.units))
print("G:\t\t\{}.".format(model.G))
print("sim t:\t\t{} = Julian day {}".format(model.t, model.t+TimeOfEphemeris))
print("Velocity:\t{} AU/yr (Earth = 2PiR AU)\n-------------------------------------".format(math.sqrt(model.particles["Earth"].vx**2 + model.particles[labels.index("Earth")].vy**2 + model.particles["Earth"].vz**2)))

model.calculate_orbits()
print("Crosscheck Io:\nx={}, y={}, z={}, vx= {}, vy={}, vz={}".format(model.particles["Io"].x, model.particles["Io"].y, model.particles["Io"].z,model.particles["Io"].vx, model.particles["Io"].vy, model.particles["Io"].vz))
print("a={}, e={}, inc={}, P={}, Omega={}, omega={}, f={}".format(model.particles["Io"].a, model.particles["Io"].e, model.particles["Io"].inc, model.particles["Io"].P*365.25, model.particles["Io"].Omega, model.particles["Io"].omega, model.particles["Io"].f))
print('Save Initial conditions as {}...'.format(BinFile))
model.save(BinFile)

print('Check...\nDelete model...')
del model
print('Load model ...')
model1 = rb.Simulation.from_file(BinFile)
model1.status()
print('Units:\t\t{}'.format(model1.units))
print("sim t:\t\t{} = Julian day {}".format(model1.t, model1.t+TimeOfEphemeris))
print("G:\t\t\t{}.".format(model1.G))
print("Velocity:\t{} AU/yr (Earth = 2PiR AU)\n-------------------------------------".format(math.sqrt(model1.particles["Earth"].vx**2 + model1.particles[labels.index("Earth")].vy**2 + model1.particles["Earth"].vz**2)))

model1.calculate_orbits()
print("Crosscheck Io:\nx={}, y={}, z={}, vx= {}, vy={}, vz={}".format(model1.particles["Io"].x, model1.particles["Io"].y, model1.particles["Io"].z,model1.particles["Io"].vx, model1.particles["Io"].vy, model1.particles["Io"].vz))
print("a={}, e={}, inc={}, P={}, Omega={}, omega={}, f={}".format(model1.particles["Io"].a, model1.particles["Io"].e, model1.particles["Io"].inc, model1.particles["Io"].P*365.25, model1.particles["Io"].Omega, model1.particles["Io"].omega, model1.particles["Io"].f))


