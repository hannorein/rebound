# Import the rebound module
import rebound

# Add particle to rebound
rebound.add( m=1. )                
rebound.add( m=1e-3, a=1., e=0.1 ) # Planet 1
rebound.add( a=1.4, e=0.1 )       # Massless test particle

# Output orbits in Jacobi coordinates
for o in rebound.calculate_orbits(): print(o)

# Output orbits in Heliocentric coordinates
for o in rebound.calculate_orbits(heliocentric=True): print(o)

# Output cartesian coordinates
for p in rebound.particles: 
    print(p)
