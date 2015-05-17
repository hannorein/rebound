# Import the rebound module
import rebound

# Add particles
# We work in units where G=1.  
rebound.add(m=1. )                  # Test particle
rebound.add(m=1e-3,x=1.,vy=1. )     # Planet

# Move particles so that the center of mass is (and stays) at the origin  
rebound.move_to_com()

# You can provide a function, written in python to REBOUND.
# This function gets called every time the forces are evaluated.
# Simple add any any additional (non-gravitational) forces to the 
# particle accelerations. Here, we add a simple drag force. This 
# will make the planet spiral into the star.
ps = rebound.particles
def dragforce():
    dragcoefficient = 1e-2
    for p in ps:
        p.ax += -dragcoefficient * p.vx
        p.ay += -dragcoefficient * p.vy
        p.az += -dragcoefficient * p.vz

# Tell rebound which function to call
rebound.additional_forces = dragforce

# Integrate until t=100 (roughly 16 orbits at 1 AU) 
rebound.integrate(100.)

# Output something at the end (the planet will be at ~0.1 AU)
for p in ps:
    print(p.x, p.y, p.z)
