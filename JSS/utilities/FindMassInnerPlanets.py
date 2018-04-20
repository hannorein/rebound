import rebound as rb
import datetime
import re

MASS_SOLAR		= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("10"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_SOLAR     /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_MERCURY  	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("199"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_MERCURY   /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_VENUS	 	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("299"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_VENUS	   /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_EARTH		= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("3"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_EARTH     /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_MARS 		= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("4"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_MARS      /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

MASS_JUPITER 	= float(re.search(r"BODY%d\_GM .* \( *([\.E\+\-0-9]+ *)\)"%int("599"), rb.horizons.HORIZONS_MASS_DATA).group(1))
MASS_JUPITER   /= rb.horizons.Gkmkgs # divide by G (horizons masses give GM. units of km^3/kg/s^2)

# set bodies for model (505 is Amalthea) and set rebound integrator conditions
labels 				= ["Sun", "Mercury", "Venus", "Earth", "Mars", "599"] 
		  
model 				= rb.Simulation()
model.units 		= ("yr", "au", "Msun")	# chosen yr instead of yr/2pill
model.integrator 	= "whfast"				# whfast integrator (no collision)

model.add(labels, date='2018-03-30 12:00') 	# Julian day 2458208

# Set particle hash to address particles by their name rather than their index

for i in range(0,len(model.particles),1):
	model.particles[i].hash=labels[i]

print('\nMass Sun:\t{0:.30f}\nMass Mercury:\t{1:.30f}\nMass Venus:\t{2:.30f}\nMass Earth:\t{3:.30f}\nMass Mars:\t{4:.30f}\nMass Jupiter:\t{5: .30f}'.format(model.particles['Sun'].m,model.particles['Mercury'].m,model.particles['Venus'].m,model.particles['Earth'].m,model.particles['Mars'].m,model.particles['599'].m))
print('---------------------------------------------------------\nTotal:\t\t{0:.30f}'.format(model.particles['Sun'].m+model.particles['Mercury'].m+model.particles['Venus'].m+model.particles['Earth'].m+model.particles['Mars'].m,))

print('\nMass Sun:\t{0:.30f}\nMass Mercury:\t{1:.30f}\nMass Venus:\t{2:.30f}\nMass Earth:\t{3:.30f}\nMass Mars:\t{4:.30f}\nMass Jupiter:\t{5: .30f}'.format(MASS_SOLAR, MASS_MERCURY, MASS_VENUS, MASS_EARTH, MASS_MARS,MASS_JUPITER))
print('---------------------------------------------------------\nTotal:\t\t{0:.30f}'.format((MASS_SOLAR+MASS_MERCURY+MASS_VENUS+MASS_EARTH+MASS_MARS)/MASS_SOLAR))

print('\nMass Sun:\t{0:.30f}\nMass Mercury:\t{1:.30f}\nMass Venus:\t{2:.30f}\nMass Earth:\t{3:.30f}\nMass Mars:\t{4:.30f}\nMass Jupiter:\t{5: .30f}'.format((MASS_SOLAR/MASS_SOLAR)-model.particles['Sun'].m,

(MASS_MERCURY/MASS_SOLAR)-model.particles['Mercury'].m, 
(MASS_VENUS/MASS_SOLAR)-model.particles['Venus'].m, 
(MASS_EARTH/MASS_SOLAR)-model.particles['Earth'].m, 
(MASS_MARS/MASS_SOLAR)-model.particles['Mars'].m, 
(MASS_JUPITER/MASS_SOLAR)-model.particles['599'].m))
print('---------------------------------------------------------\nTotal:\t\t{0:.30f}'.format((MASS_SOLAR+MASS_MERCURY+MASS_VENUS+MASS_EARTH+MASS_MARS)/MASS_SOLAR))
