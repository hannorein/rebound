# Import the rebound module
import rebound

# Create a REBOUND simulation
sim = rebound.Simulation()

# Set variables (defaults are G=1, t=0, dt=0.01)
k = 0.01720209895       # Gaussian constant 
sim.G = k*k         # Gravitational constant

# Setup particles (data taken from NASA Horizons)
# This could also be easily read in from a file.
sim.add( m=1.00000597682, x=-4.06428567034226e-3, y=-6.08813756435987e-3, z=-1.66162304225834e-6, vx=+6.69048890636161e-6, vy=-6.33922479583593e-6, vz=-3.13202145590767e-9)   # Sun
sim.add( m=1./1047.355,   x=+3.40546614227466e+0, y=+3.62978190075864e+0, z=+3.42386261766577e-2, vx=-5.59797969310664e-3, vy=+5.51815399480116e-3, vz=-2.66711392865591e-6)   # Jupiter
sim.add( m=1./3501.6,     x=+6.60801554403466e+0, y=+6.38084674585064e+0, z=-1.36145963724542e-1, vx=-4.17354020307064e-3, vy=+3.99723751748116e-3, vz=+1.67206320571441e-5)   # Saturn
sim.add( m=1./22869.,     x=+1.11636331405597e+1, y=+1.60373479057256e+1, z=+3.61783279369958e-1, vx=-3.25884806151064e-3, vy=+2.06438412905916e-3, vz=-2.17699042180559e-5)   # Uranus
sim.add( m=1./19314.,     x=-3.01777243405203e+1, y=+1.91155314998064e+0, z=-1.53887595621042e-1, vx=-2.17471785045538e-4, vy=-3.11361111025884e-3, vz=+3.58344705491441e-5)   # Neptune
sim.add( m=0,             x=-2.13858977531573e+1, y=+3.20719104739886e+1, z=+2.49245689556096e+0, vx=-1.76936577252484e-3, vy=-2.06720938381724e-3, vz=+6.58091931493844e-4)   # Pluto

# Set the center of momentum to be at the origin
sim.move_to_com()

# timestep counter
steps = 0 
# Integrate until t=1e4 (unit of time in this example is days)
for i in range(10):
    t = 1.0e3*i
    sim.integrate(t)
    # Print particle positions 
    for p in sim.particles:
        #     time       x    y    z 
        print(sim.t, p.x, p.y, p.z)

