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
        self.sim.particles[2].hash = rebound.hash("earth")
        self.sim.particles[3].hash = rebound.hash("jupiter")

    def tearDown(self):
        self.sim = None

    def test_add_hash(self):
        self.assertEqual(self.sim.particles[0].hash, 0)
        self.assertEqual(self.sim.particles[2].hash, 1424801690)

    def test_get_particle_by_string(self):
        self.assertAlmostEqual(self.sim.get_particle_by_string("earth").a, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle_by_string("jupiter").a, 5., delta=1e-15)
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.get_particle_by_string("venus")

    def test_removing_particles(self):
        self.sim.remove(index=1)
        self.assertAlmostEqual(self.sim.get_particle_by_string("earth").a, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle_by_string("jupiter").a, 5., delta=1e-15)

    def test_adding_particles(self):
        self.assertAlmostEqual(self.sim.get_particle_by_string("earth").a, 1., delta=1e-15)
        self.assertEqual(self.sim.N_particle_lookup_table, 2)
        self.sim.add(a=10.)
        self.sim.particles[4].hash = rebound.hash("saturn")
        self.assertAlmostEqual(self.sim.get_particle_by_string("jupiter").a, 5., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle_by_string("saturn").a, 10., delta=1e-15)
        self.assertAlmostEqual(self.sim.get_particle_by_string("earth").a, 1., delta=1e-15)
        self.assertEqual(self.sim.N_particle_lookup_table, 3)

if __name__ == "__main__":
    unittest.main()
