import rebound
import unittest
import os

class TestUnits(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_units_no_particle(self):
        default_units = {'length':None, 'mass':None, 'time':None}
        self.assertEqual(self.sim.units, default_units)
        new_units = {'length':'au', 'mass':'msun', 'time': 'yr2pi'}
        self.sim.units = ("au", "msun", "yr2pi")
        self.assertEqual(self.sim.units, new_units)
        
        self.sim.add(m=1.)
        with self.assertRaises(AttributeError):
            self.sim.units = ("au", "kg", "yr2pi")
    
    def test_wrong_units(self):
        with self.assertRaises(Exception):
            self.sim.units = ("au", "yr2pi")
        with self.assertRaises(Exception):
            self.sim.units = ("au", "bogusunit", "yr2pi")
    
    def test_units_with_particle(self):
        self.sim.units = ("au", "msun", "yr2pi")
        self.sim.add(m=1.)
        self.sim.add(m=1., x=1., vx=1.)
        self.sim.convert_particle_units("au", "kg", "yr2pi")
        self.assertAlmostEqual(self.sim.particles[0].m,1.988499251452493e+30, delta=1e-15)
        self.sim.convert_particle_units("m", "kg", "yr2pi")
        self.assertAlmostEqual(self.sim.particles[1].x,149597870700.0, delta=1e-15)
        self.sim.convert_particle_units("m", "kg", "s")
        self.assertAlmostEqual(self.sim.particles[1].vx,29784.691834383168, delta=1e-15)

