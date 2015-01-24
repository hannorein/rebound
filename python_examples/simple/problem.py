# Import the rebound module
import sys; sys.path.append('../')
import rebound

# Set variables (defaults are G=1, t=0, dt=0.01)
#rebound.set_G(1.)  
#rebound.set_t(0.)  
#rebound.set_dt(0.01)  

# Add particles
# All parameters omitted are set to 0 by default.
rebound.particle_add( m=1. )                  # Star
rebound.particle_add( m=1e-3, x=1., vy=1. )   # Planet

# Move particles so that the center of mass is (and stays) at the origin  
rebound.move_to_center_of_momentum()

# Integrate until t=100 (roughly 16 orbits) 
rebound.integrate(100.)

# Modify particles 
# As an example, we are reverting the velocities 
particles = rebound.particles_get()
for i in range(rebound.get_N()):
    particles[i].vx *= -1.
    particles[i].vy *= -1.
    particles[i].vz *= -1.

# Integrate another 100 time units, until t=200
rebound.integrate(200.)

# Get particles back and print positions 
# Since we're integrating forward and then backward we end up with the 
# particles exactly where they started out from (note that we moved to the
# center of momentum frame)
for i in range(rebound.get_N()):
    print(particles[i].x, particles[i].y, particles[i].z)
