import rebound
import math
import random
import matplotlib.pyplot as plt
import numpy as np

sim = rebound.Simulation()
OMEGA = 0.0001485     # [1/s]

surface_density = 250.    # kg/m^2
particle_density = 400.   # kg/m^3
sim.G = 6.67428e-11       # N m^2 / kg^2
sim.dt = 1e-3*2.*math.pi/OMEGA
sim.softening = 0.1       # [m]
boxsize = 100.            # [m]
sim.configure_box(boxsize, 2, 2, 1)
sim.N_ghost_x = 2
sim.N_ghost_y = 2
sim.N_ghost_z = 0
sim.integrator = "sei"
sim.boundary   = 3
sim.gravity    = "tree"
sim.collision  = "tree"
sim.collision_resolve = "hardsphere"
sim.ri_sei.OMEGA = OMEGA
sim.ri_sei.Q_NL = 0.0153
Q_NL = 0.1
def cor_bridges(r, v):
    eps = 0.32*pow(abs(v)*100.,-0.234)
    if eps>1.:
        eps=1.
    if eps<0.:
        eps=0.
    return eps
sim.coefficient_of_restitution = cor_bridges
def powerlaw(slope, min_v, max_v):
    y = random.random()
    pow_max = pow(max_v, slope+1.)
    pow_min = pow(min_v, slope+1.)
    return pow((pow_max-pow_min)*y + pow_min, 1./(slope+1.))

total_mass = surface_density*sim.boxsize.x*sim.boxsize.y
mass = 0.0
while mass < total_mass:
    radius = powerlaw(slope=-3, min_v=1, max_v=4)  # [m]    
    m = particle_density*4./3.*math.pi*(radius**3)
    x = random.uniform(-boxsize/2., boxsize/2.)
    y = random.uniform(-boxsize/2., boxsize/2.)
    sim.add(
        m=mass,
        r=radius,
        x=x-Q_NL*x*math.cos(OMEGA*sim.t),
        y=y-1.5*x*OMEGA*sim.t+2.0*Q_NL*x*math.sin(OMEGA*sim.t),
        z=0.0,
        vx = Q_NL*x*OMEGA*math.sin(OMEGA*sim.t),
        vy = -1.5*x*OMEGA+2.0*Q_NL*x*OMEGA*math.cos(OMEGA*sim.t), 
        vz = 0.)
    mass += m

print(len(sim.particles))

period = 2*math.pi/OMEGA
epsilon = 1e-4
x = []
vx = []
for i in range(100):
    sim.integrate(i*period/100)
    p = sim.particles[0]
    if p.y <= epsilon:
        x.append(p.x)
        vx.append(p.vx)

print(x)
print(vx)
_x = np.array(x)
_vx = np.array(vx)
plt.scatter(_x, _vx)
plt.show()
