from rebound import *

particles = []
particles.append( Particle(m=1.) )
particles.append( Particle(m=1e-3,x=1.,vy=1.) )

set_particles(particles)

while get_t() < 100.:
    step()

particles = get_particles()
for p in particles:
    print p.x, p.y, p.y

