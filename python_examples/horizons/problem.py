import rebound
# Need to specify mass as HORIZONS doesn't return physical properties in easily readable format.
# Optionally set date with argument date="2015-05-08 12:00" 
rebound.add("Sun",m=1.)
rebound.add("EMB") # Earth Moon Barycenter (sometimes stupid name convention)
rebound.add("Jupiter Barycenter",m=0.0009547919) 
rebound.add("C/2014 Q2") # Comets!
rebound.status()
rebound.set_integrator("whfast")
rebound.integrate(100.)
rebound.status()
