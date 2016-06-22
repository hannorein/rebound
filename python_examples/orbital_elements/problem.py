import rebound

# Create a rebound simulation
sim = rebound.Simulation()

# Add a particle at the origin with mass 1
sim.add(m=1.)  

# Add a particle with mass 1e-3 on a Keplerian 
# orbit around the center of mass (including all 
# particles added so far) with semi-major axis 1
sim.add(m=1e-3, a=1.)

# Add a test particle (mass=0) on a Keplerian orbit 
# around the center of mass (both particles added above)
# with a semi-major axis of 2 and eccentricity of 0.1.
# This corresponds to Jacobi coordinates.
sim.add(a=1., e=0.1)

# Move all particles to the-center-of-momentum frame.
sim.move_to_com()

# Print the resulting cartesian coordinates.
for p in sim.particles:
    print(p.m, p.x, p.y, p.z, p.vx, p.vy, p.vz)

# Integrate for 100 time units
sim.integrate(100.)

