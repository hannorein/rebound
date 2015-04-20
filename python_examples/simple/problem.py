# Import the rebound module
import rebound

def print_orbits(): 
    print("time: %f"%rebound.get_t())
    # Get particles back and print orbital parameters 
    # Since we're integrating forward and then backward we end up with the 
    # particles exactly where they started out from 
    for i in range(1,rebound.get_N()):
        print(particles[i].get_orbit())

# Set variables (defaults are G=1, t=0, dt=0.01)
#rebound.set_G(1.)  
#rebound.set_t(0.)  
#rebound.set_dt(0.5)  

# Choose integrator (default is IAS15)
#rebound.set_integrator("mikkola")

# Set integrator options
#rebound.set_integrator_mikkola_corrector(7) # This enables a 7th order symplectic corrector.

# Add particles
# This example shows two different ways to add particles.
# All parameters omitted are set to 0 by default.
sun = rebound.Particle( m=1. )                          # Star
rebound.add_particle( sun )                
rebound.add_particle( primary=sun, m=1e-3, a=1, e=0.1 ) # Planet 1
rebound.add_particle( primary=sun, a=1.4, e=0.1 )       # Massless test particle

# Move particles so that the center of mass is (and stays) at the origin  
rebound.move_to_center_of_momentum()

# Get a pointer to the particle structure
particles = rebound.get_particles()

# Print initial conditions
print_orbits()

# Integrate until t=100 (1 orbit=2pi in these units, thus roughly 16 orbits) 
rebound.integrate(100.)

# Modify particles 
# As an example, we are reverting the velocities of all particles 
for i in range(rebound.get_N()):
    particles[i].vx *= -1.
    particles[i].vy *= -1.
    particles[i].vz *= -1.

# Print new particle orbits
print_orbits()

# Integrate another 100 time units, until t=200
rebound.integrate(200.)

# Print final orbital parameters 
# Since we're integrating forward and then backward we end up with the 
# particles exactly where they started out from 
print_orbits()
