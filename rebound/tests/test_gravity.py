import rebound
import unittest
import math
import numpy as np

class TestGravity(unittest.TestCase):
    
    def test_testparticle_0(self):
        sim = rebound.Simulation()
        sim.testparticle_type = 0;
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        sim.integrate(10.)
        x0 = sim.particles[0].x
        self.assertEqual(x0, 0.)

    def test_testparticle_1(self):
        sim = rebound.Simulation()
        sim.testparticle_type = 1;
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        sim.integrate(10.)
        x0 = sim.particles[0].x
        self.assertNotEqual(x0, 0.)
    
    def test_testparticle_comp_0(self):
        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.testparticle_type = 0;
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        sim.integrate(10.)
        x0 = sim.particles[0].x
        self.assertEqual(x0, 0.)

    def test_testparticle_comp_1(self):
        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.testparticle_type = 1;
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        sim.integrate(10.)
        x0 = sim.particles[0].x
        self.assertNotEqual(x0, 0.)


if __name__ == "__main__":
    unittest.main()
