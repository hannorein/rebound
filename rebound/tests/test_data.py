import rebound
import unittest

class TestData(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_add_solar_system(self):
        self.sim.add("solar system")
        self.assertEqual(self.sim.N,9)
        self.assertAlmostEqual(self.sim.particles[8].x,29.39189844361884,delta=1e-15)
        self.assertAlmostEqual(self.sim.particles[8].a,30.076046335613285,delta=1e-15)
    
    def test_add_outer_solar_system(self):
        self.sim.add("outer solar system")
        self.assertEqual(self.sim.N,5)
        self.assertAlmostEqual(self.sim.particles[4].x,29.39189844361884,delta=1e-15)
        self.assertAlmostEqual(self.sim.particles[4].a,30.077771073224017,delta=1e-15)

    
if __name__ == "__main__":
    unittest.main()
