import rebound
import unittest
from ctypes import c_uint32

class TestHash(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(a=0.1)
        self.sim.add(a=1.)
        # Adding particles with names in different ways
        # to test memory management
        self.sim.add(a=5.,name="jupiter")  
        self.sim.particles[2].name = "earth"
        p = rebound.Particle(x=1,vy=1,name="💫")
        self.sim.add(p)

    def tearDown(self):
        self.sim = None

    def test_compare(self):
        self.assertEqual(self.sim.particles[0].name, None)
        self.assertEqual(self.sim.particles[3].name, "jupiter")

    def test_find_particles_by_name(self):
        self.assertAlmostEqual(self.sim.particles["earth"].a, 1., delta=1e-15)
        self.assertAlmostEqual(self.sim.particles["jupiter"].a, 5., delta=1e-15)
        self.assertEqual(self.sim.particles["💫"].m, 0.0)
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["venus"]

    def test_removing_particles(self):
        self.sim.add(a=30.)
        self.sim.remove(0)
        self.sim.remove(0)
        self.sim.remove("earth", keep_sorted=0)
        self.sim.remove("💫", keep_sorted=1)
        self.assertEqual(self.sim.particles["jupiter"].name, "jupiter")
        self.assertEqual(self.sim.N, 2)
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["earth"]
        self.sim.remove("jupiter")
        with self.assertRaises(rebound.ParticleNotFound):
            self.sim.particles["jupiter"]

if __name__ == "__main__":
    unittest.main()
