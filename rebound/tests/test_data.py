import rebound
import unittest
import os
import rebound.data as data

class TestData(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_add_outer_solar_system(self):
        data.add_outer_solar_system(self.sim)
        self.assertEqual(self.sim.N,6)
        self.assertAlmostEqual(self.sim.particles[5].x,-21.3858977531573,delta=1e-15)
        self.assertAlmostEqual(self.sim.particles[5].a,39.485206935092286,delta=1e-15)

