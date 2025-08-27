import rebound
import unittest
import sys
import warnings
import math

class TestOrbitHierarchy(unittest.TestCase):
    
    def test_solar_system(self):
        sim = rebound.Simulation()
        sim.add("solar system")
        oh = sim.orbit_hierarchy()
        self.assertEqual(True, oh.is_jacobi())
        self.assertEqual(True, oh.is_jacobi_ordered(sim))
        
        sim = rebound.Simulation()
        sim.add("outer solar system")
        oh = sim.orbit_hierarchy()
        self.assertEqual(True, oh.is_jacobi())
        self.assertEqual(True, oh.is_jacobi_ordered(sim))


    def test_not_jacobi(self):
        sim = rebound.Simulation()
        sim.add("solar system")
        sim.add(m=1e-8, a=0.001, primary=sim.particles[3]) # moon
        oh = sim.orbit_hierarchy()
        self.assertEqual(False, oh.is_jacobi())
        self.assertEqual(False, oh.is_jacobi_ordered(sim))
    
    def test_out_of_order(self):
        sim = rebound.Simulation()
        sim.add("solar system")
        sim.add(m=1e-8, a=0.001,primary=sim.particles[0]) # close in planet
        oh = sim.orbit_hierarchy()
        self.assertEqual(True, oh.is_jacobi())
        self.assertEqual(False, oh.is_jacobi_ordered(sim))
        
if __name__ == "__main__":
    unittest.main()

