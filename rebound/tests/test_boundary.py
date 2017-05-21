import rebound
import unittest

class TestBoundary(unittest.TestCase):
    
    def test_open(self):
        sim = rebound.Simulation()
        sim.boundary = "open"
        sim.configure_box(10.)
        sim.add(m=0.1,x=1., vx=5.0)
        sim.add(m=0.1,x=-1., vx=-5.0)
        sim.add(m=0.1,y=1., vx=6.0)
        sim.add(m=0.1,x=-1., y=-1., vx=-3., vy=-3.)
        self.assertEqual(sim.N,4)
        sim.integrate(1.)
        self.assertEqual(sim.N,1)
        with self.assertRaises(rebound.NoParticles):
            sim.integrate(2.)
    
    def test_periodic(self):
        sim = rebound.Simulation()
        sim.boundary = "periodic"
        sim.configure_box(10.)
        sim.add(m=0.1,x=1., vx=5.0, vy=15.1, vz=26.)
        sim.integrate(1.)
        self.assertAlmostEqual(sim.particles[0].x,-4,delta=1e-16)
        sim.integrate(2.)
        self.assertAlmostEqual(sim.particles[0].x,1,delta=1e-16)
        self.assertEqual(sim.N,1)
    
    
if __name__ == "__main__":
    unittest.main()
