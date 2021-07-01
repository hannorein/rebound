import rebound
import unittest
import math
import rebound.data
import warnings
import copy

class TestIntegratorTES(unittest.TestCase):
    def test_particle_pass_through(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3,a=1.12313)
        sim.add(m=1e-3,a=2.32323)
        sim.move_to_com()
        sim.dt = 0.25
        sim.integrator = "tes"
        pos = sim.particles[1].x
        sim.integrate(1)

        self.assertEqual(pos, sim.particles[1].x)




if __name__ == "__main__":
    unittest.main()
    # sim = rebound.Simulation()
    # sim.add(m=1.)
    # sim.add(m=1e-3,a=1.12313)
    # sim.add(m=1e-3,a=2.32323)
    # sim.move_to_com()
    # sim.dt = 0.5
    # sim.integrator = "tes"
    # e0 = sim.calculate_energy()
    # sim.integrate(1)
    # e1 = sim.calculate_energy()
    # self.assertLess(math.fabs((e0-e1)/e1),eps)
    