import rebound
import unittest
import math
import random

class TestShearingSheet(unittest.TestCase):
    
    def test_saturnsrings(self):
        sim = rebound.Simulation()
        OMEGA = 0.00013143527     # [1/s]
        sim.ri_sei.OMEGA = OMEGA
        surface_density = 400.    # kg/m^2
        particle_density = 400.   # kg/m^3
        sim.G = 6.67428e-11       # N m^2 / kg^2
        sim.dt = 1e-3*2.*math.pi/OMEGA
        sim.softening = 0.2       # [m]
        boxsize = 50.            # [m]
        sim.configure_box(boxsize)
        sim.N_ghost_x = 2
        sim.N_ghost_y = 2
        sim.integrator = "sei"
        sim.boundary   = "shear"
        sim.gravity    = "tree"
        sim.collision  = "tree"
        sim.collision_resolve = "hardsphere"
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
        total_mass = 0.
        while total_mass < surface_density*(boxsize**2):
            radius = powerlaw(slope=-3, min_v=1, max_v=4)  # [m]    
            mass = particle_density*4./3.*math.pi*(radius**3)
            x = random.uniform(-boxsize/2., boxsize/2.)
            sim.add(
                m=mass,
                r=radius,
                x=x,
                y=random.uniform(-boxsize/2., boxsize/2.),
                z=random.normalvariate(mu=0.0, sigma=1.0),
                vx = 0.,
                vy = -3./2.*x*OMEGA, 
                vz = 0.)
            total_mass += mass
        self.assertGreater(sim.N,50)
        sim.integrate(2.*math.pi/OMEGA)
        self.assertGreater(sim.collisions_log_n,1000)
        Nbefore = sim.N
        sim.remove(0,keep_sorted=0)
        sim.update_tree()
        self.assertEqual(Nbefore-1,sim.N)
        with self.assertRaises(RuntimeError):
            sim.remove(0,keep_sorted=1)
        self.assertNotEqual(sim.ri_sei._lastdt,0.0)
        sim.reset_integrator()
        self.assertEqual(sim.ri_sei._lastdt,0.0)

if __name__ == "__main__":
    unittest.main()
