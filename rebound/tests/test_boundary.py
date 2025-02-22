import rebound
import unittest

class TestBoundary(unittest.TestCase):
    
    def test_open(self):
        for gravity in ["basic", "tree"]:
            sim = rebound.Simulation()
            sim.boundary = "open"
            sim.gravity  = gravity
            sim.configure_box(10.)
            sim.add(m=0.1,x=1., vx=5.0)
            sim.add(m=0.1,x=-1., vx=-5.0)
            sim.add(m=0.1,y=1., vx=6.0)
            sim.add(m=0.1,x=-1., y=-1., vx=-3., vy=-3.)
            sim.add(m=0.1,z=1., vz=5.0)
            sim.add(m=0.1,z=-1., vz=-5.0)
            self.assertEqual(sim.N,6)
            sim.integrate(1.)
            self.assertEqual(sim.N,1)
            with self.assertRaises(rebound.NoParticles):
                sim.integrate(2.)
    
    def test_periodic(self):
        sim = rebound.Simulation()
        sim.boundary = "periodic"
        sim.configure_box(10.)
        sim.add(m=0.0,x=1., vx=5.0, vy=15.1, vz=26.)
        sim.add(m=0.0,x=-1., vx=-5.0, vy=-15.1, vz=-26.)
        sim.integrate(1.)
        self.assertAlmostEqual(sim.particles[0].x,-4,delta=1e-16)
        self.assertAlmostEqual(sim.particles[0].y,-4.9,delta=1e-14)
        self.assertAlmostEqual(sim.particles[0].z,-4.0,delta=1e-16)
        self.assertAlmostEqual(sim.particles[1].x,4,delta=1e-16)
        sim.integrate(2.)
        self.assertAlmostEqual(sim.particles[0].x,1,delta=1e-16)
        self.assertAlmostEqual(sim.particles[0].y,0.2,delta=1e-14)
        self.assertAlmostEqual(sim.particles[0].z,2.0,delta=1e-16)
        self.assertAlmostEqual(sim.particles[1].x,-1,delta=1e-16)
        self.assertEqual(sim.N,2)
    
    def test_shear_vertical(self):
        # Note: not physical
        sim = rebound.Simulation()
        sim.ri_sei.OMEGA = 1
        sim.configure_box(1)
        sim.N_ghost_x = 2
        sim.N_ghost_y = 2
        sim.integrator = "sei"
        sim.boundary   = "shear"
        sim.add(z=0.45,vz=100.0) 
        sim.add(z=-0.45,vz=-100.0) 
        sim.integrate(1)
        self.assertEqual(2,sim.N)

    def test_add_outside(self):
        for boundary in ["open", "shear", "periodic"]:
            sim = rebound.Simulation()
            sim.boundary = boundary
            sim.configure_box(10.)
            with self.assertRaises(RuntimeError):
                sim.add(x=6.)
            with self.assertRaises(RuntimeError):
                sim.add(x=-6.)
            with self.assertRaises(RuntimeError):
                sim.add(y=6.)
            with self.assertRaises(RuntimeError):
                sim.add(y=-6.)
            with self.assertRaises(RuntimeError):
                sim.add(z=6.)
            with self.assertRaises(RuntimeError):
                sim.add(z=-6.)
    
    
if __name__ == "__main__":
    unittest.main()
