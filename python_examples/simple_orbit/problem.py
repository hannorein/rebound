# Import the rebound module
import rebound

# Create Simulation object
sim = rebound.Simulation()
# Add particle to rebound
sim.add( m=1. )                
sim.add( m=1e-3, a=1., e=0.1 ) # Planet 1
sim.add( a=1.4, e=0.1 )       # Massless test particle

# Output orbits in Jacobi coordinates
for o in sim.orbits(): print(o)

# Output orbits in Heliocentric coordinates
for o in sim.orbits(primary=sim.particles[0]): print(o)

# Output cartesian coordinates
for p in sim.particles: 
    print(p)
