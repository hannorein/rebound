# Import the rebound module
import rebound

# Add particle to rebound
rebound.add( m=1. )                
rebound.add( m=1e-3, a=1., e=0.1 ) # Planet 1
rebound.add( a=1.4, e=0.1 )       # Massless test particle

# Move particles so that the center of mass is (and stays) at the origin  
rebound.move_to_center_of_momentum()

# Integrate until t=100 (1 orbit=2pi in these units, thus roughly 16 orbits) 
#rebound.integrate(100.)

for o in rebound.get_orbits():
    print o
