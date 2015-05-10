import rebound
rebound.add("Sun")
#rebound.add("Mercury")
#rebound.add("Earth") 
#rebound.add("C/2014 Q2") # Comets!
rebound.status()
rebound.set_integrator("whfast")
rebound.set_dt(0.01)
rebound.integrate(100.)
rebound.status()
