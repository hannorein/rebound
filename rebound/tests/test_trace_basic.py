import rebound
import unittest
import numpy as np


class TestIntegratorTraceHarmonic(unittest.TestCase):
   
    def test_trace_harmonic_with_nbody(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1);
        sim.add(m=1e-3,a=1.1);
        sim.integrator = "TRACE"

        sim.step()
        
        self.assertEqual(sim.ri_trace._current_C,0)
        self.assertEqual(sim.ri_trace._current_Ks[1],0)


if __name__ == "__main__":
    unittest.main()
