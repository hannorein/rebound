import rebound
import unittest
import datetime

class TestHorizons(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_earth(self):
        self.sim.add("Earth",date="2000-01-01 00:00")
        self.assertAlmostEqual(self.sim.particles[0].x,-0.16855044899082616,delta=1e-10)
        self.assertAlmostEqual(self.sim.particles[0].m,3.0404326480226416e-06,delta=1e-15)
    
    def test_notfound(self):
        with self.assertRaises(Exception):
            self.sim.add("BogusPlanet",date="2000-01-01 00:00")


