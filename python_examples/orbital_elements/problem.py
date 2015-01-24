import sys; sys.path.append('../')
import rebound

# Add a particle at the origin with mass 1
rebound.particle_add(m=1.)  

# Add a particle with mass 1e-3 on a Keplerian 
# orbit around the center of mass (including all 
# particles added so far) with semi-major axis 1
rebound.particle_add(m=1e-3, a=1.)

# Add a test particle (mass=0) on a Keplerian orbit 
# around the center of mass (both particles added above)
# with a semi-major axis of 2 and eccentricity of 0.1.
# This corresponds to Jacobi coordinates.
rebound.particle_add(a=1., e=0.1)

# Move all particles to the-center-of-momentum frame.
rebound.move_to_center_of_momentum()

# Print the resulting cartesian coordinates.
p = rebound.particles_get()
for i in range(rebound.get_N()):
    print(p[i].m, p[i].x, p[i].y, p[i].z, p[i].vx, p[i].vy, p[i].vz)

# Integrate for 100 time units
rebound.integrate(100.)

