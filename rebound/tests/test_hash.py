import rebound
import unittest
import math 
from ctypes import c_uint32

class TestHash(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=0.1)
        self.sim.add(a=1.)
        self.sim.add(a=5.)
        self.sim.particles[2].hash = "earth"
        self.sim.particles[3].hash = "jupiter"
        self.uint32_max = 4294967295

    def tearDown(self):
        self.sim = None

    def test_rebound_hash(self):
        self.assertNotEqual(rebound.hash("foo").value, rebound.hash("bar").value)
        self.assertEqual(rebound.hash("earth").value, 1424801690)
        
    def test_add_hash(self):
        self.assertEqual(self.sim.particles[0].hash.value, self.uint32_max)
        self.assertEqual(self.sim.particles[2].hash.value, 1424801690)

    def test_particles(self):
        self.assertAlmostEqual(self.sim.particles["earth"].a, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.particles["jupiter"].a, 5., delta=1e-15)
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["venus"]

    def test_removing_particles(self):
        self.sim.add(a=30.)
        self.sim.remove(0)
        self.sim.remove(0)
        self.sim.remove(hash=rebound.hash("earth"), keepSorted=0)
        self.assertEqual(self.sim.particles["jupiter"].hash.value, rebound.hash("jupiter").value)
        self.assertEqual(self.sim.N, 2)
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["earth"]
        self.sim.remove(hash="jupiter")
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["jupiter"]

    def test_adding_particles(self):
        self.assertAlmostEqual(self.sim.particles["earth"].a, 1., delta=1e-15)
        self.assertEqual(self.sim.N_lookup, 2)
        self.sim.add(a=10.)
        self.sim.particles[4].hash = "saturn"
        self.assertAlmostEqual(self.sim.particles["jupiter"].a, 5., delta=1e-15)
        self.assertAlmostEqual(self.sim.particles["saturn"].a, 10., delta=1e-15)
        self.assertAlmostEqual(self.sim.particles["earth"].a, 1., delta=1e-15)
        self.assertEqual(self.sim.N_lookup, 3)

    def test_add(self):
        self.sim.add(a=7., hash = "planet 9")
        self.assertEqual(rebound.hash("planet 9").value, self.sim.particles["planet 9"].hash.value)

class TestEmptyLookup(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=0.1)
        self.sim.add(a=1.)
        self.sim.add(a=5.)

    def tearDown(self):
        self.sim = None

    def test_empty_lookup(self):
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["venus"]

class TestZeroHash(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()

    def tearDown(self):
        self.sim = None

    def test_zero_first_hash(self):
        self.sim.add(m=1., hash=c_uint32(0))
        self.sim.add(m=2., hash=c_uint32(1))
        self.assertAlmostEqual(self.sim.particles[c_uint32(1)].m, 2., delta=1e-15)
        self.assertEqual(self.sim.N_lookup, 2)
        self.sim.add(m=3., hash="jupiter")
        self.assertAlmostEqual(self.sim.particles[c_uint32(0)].m, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.particles["jupiter"].m, 3., delta=1e-15)
        self.assertEqual(self.sim.N_lookup, 3)

if __name__ == "__main__":
    unittest.main()
