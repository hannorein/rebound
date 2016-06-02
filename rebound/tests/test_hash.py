import rebound
import unittest
import math 

class TestHash(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=0.1)
        self.sim.add(a=1.)
        self.sim.add(a=5.)
        self.sim.particles[2].hash = self.sim.get_hash("earth")
        self.sim.particles[3].hash = self.sim.get_hash("jupiter")

    def tearDown(self):
        self.sim = None

    def test_add_hash(self):
        self.assertEqual(self.sim.particles[0].hash, 0)
        self.assertEqual(self.sim.particles[2].hash, 1424801690)

    def test_get_particle(self):
        self.assertAlmostEqual(self.sim.get_particle("earth").a, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle("jupiter").a, 5., delta=1e-15)
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.get_particle("venus")

    def test_removing_particles(self):
        self.sim.remove(index=1)
        self.assertAlmostEqual(self.sim.get_particle("earth").a, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle("jupiter").a, 5., delta=1e-15)

    def test_adding_particles(self):
        self.assertAlmostEqual(self.sim.get_particle("earth").a, 1., delta=1e-15)
        self.assertEqual(self.sim.N_lookup, 2)
        self.sim.add(a=10.)
        self.sim.particles[4].hash = self.sim.get_hash("saturn")
        self.assertAlmostEqual(self.sim.get_particle("jupiter").a, 5., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle("saturn").a, 10., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle("earth").a, 1., delta=1e-15)
        self.assertEqual(self.sim.N_lookup, 3)

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
            self.sim.get_particle("venus")

class TestAssignedHashes(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        for x in range(100):
            self.sim.add(x=x)
            self.sim.particles[-1].hash = self.sim.get_hash()

    def tearDown(self):
        self.sim = None

    def test_assigned_hashes(self):
        hash53 = self.sim.particles[53].hash
        hash54 = self.sim.get_particle(hash=hash53+1).hash
        self.assertEqual(self.sim.N_lookup, 100)
        self.assertEqual(hash54, hash53+1)

if __name__ == "__main__":
    unittest.main()
