import rebound
import unittest
import math
import rebound.data

class TestIntegrator(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(self.sim)
        self.sim.move_to_com()
    
    def tearDown(self):
        self.sim = None
    
    def test_leapfrog_order_2(self):
        self.sim.integrator = "leapfrog"
        self.sim.ri_leapfrog.order = 2
        self.sim.dt = 1
        e0 = self.sim.energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3)
        e1 = self.sim.energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-3)
    
    def test_leapfrog_order_4(self):
        self.sim.integrator = "leapfrog"
        self.sim.ri_leapfrog.order = 4
        self.sim.dt = 1
        e0 = self.sim.energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3)
        self.sim.step()
        e1 = self.sim.energy()
        self.assertLess(math.fabs((e0-e1)/e1),3e-6)

    def test_leapfrog_order_6(self):
        self.sim.integrator = "leapfrog"
        self.sim.ri_leapfrog.order = 6
        self.sim.dt = 1
        e0 = self.sim.energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3)
        self.sim.step()
        e1 = self.sim.energy()
        self.assertLess(math.fabs((e0-e1)/e1),2e-9)
    
    def test_leapfrog_order_8(self):
        self.sim.integrator = "leapfrog"
        self.sim.ri_leapfrog.order = 8
        self.sim.dt = 1
        e0 = self.sim.energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3)
        self.sim.step()
        e1 = self.sim.energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e13)
    
    def test_leapfrog_wrong_order(self):
        self.sim.integrator = "leapfrog"
        self.sim.ri_leapfrog.order = 42
        with self.assertRaises(RuntimeError):
            self.sim.step()

if __name__ == "__main__":
    unittest.main()
