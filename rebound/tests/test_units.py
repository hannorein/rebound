import rebound
import unittest
import os

class TestSimulationConstructor(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_units(self):
        default_units = {'length':None, 'mass':None, 'time':None}
        self.assertEqual(self.sim.units, default_units)
        new_units = {'length':'au', 'mass':'msun', 'time': 'yr2pi'}
        self.sim.units = ("au", "msun", "yr2pi")
        self.assertEqual(self.sim.units, new_units)
        
        self.sim.add(m=1.)
        with self.assertRaises(AttributeError):
            self.sim.units = ("au", "kg", "yr2pi")

