# Import the rebound module
import sys; sys.path.append('../')
import rebound
from rebound import Particle

# Add particles
# We work in units where G=1.  
rebound.particle_add( Particle(m=1.) )                  # Test particle
rebound.particle_add( Particle(m=1e-3,x=1.,vy=1.) )     # Planet

# Move particles so that the center of mass is (and stays) at the origin  
rebound.move_to_center_of_momentum()

# You can provide a function, written in python to REBOUND.
# This function gets called every time the forces are evaluated.
# Simple add any any additional (non-gravitational) forces to the 
# particle accelerations. Here, we add a simple drag force. This 
# will make the planet spiral into the star.
particles = rebound.particles_get()                     # Pointer to the particle structure
N = rebound.get_N()                 
def dragforce():
    dragcoefficient = 1e-2
    for i in range(N):
        particles[i].ax += -dragcoefficient * particles[i].vx
        particles[i].ay += -dragcoefficient * particles[i].vy
        particles[i].az += -dragcoefficient * particles[i].vz

# Tell rebound which function to call
rebound.set_additional_forces(dragforce)

# Integrate until t=100 (roughly 16 orbits at 1 AU) 
rebound.integrate(100.)

# Output something at the end (the planet will be at ~0.1 AU)
for i in range(N):
    print(particles[i].x, particles[i].y, particles[i].z)
