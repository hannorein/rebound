import sys
sys.path.append('../')
from rebound import *

particles_add( Particle(m=1.) )
particles_add( Particle(m=1e-3,x=1.,vy=1.) )

move_to_center_of_momentum()

p0 = get_particle(0)
p1 = get_particle(1)
print p1.x, p1.y, p1.z
print p1.x-p0.x, p1.y-p0.y, p1.z-p0.z

integrate(100.)

particles = get_particles()
p0 = particles[0]
p1 = particles[1]
print p1.x, p1.y, p1.z
print p1.x-p0.x, p1.y-p0.y, p1.z-p0.z

