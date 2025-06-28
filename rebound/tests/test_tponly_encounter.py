import rebound
import unittest
import sys
import warnings
import math
import random
from datetime import datetime

class TestTPOnlyEncounter(unittest.TestCase):
    
    def test_tponly_encounter0_mercurius(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-2,a=1)
        sim.add(m=1e-2,a=2)
        sim.move_to_com()
        sim.N_active = sim.N
        sim.N_active = sim.N
        sim.integrator = "mercurius"
        sim.ri_mercurius.tponly_encounter = 0 # Default
        sim.dt = 0.1
        sim2 = sim.copy()
        
        random.seed(10)
        for i in range(10):
            sim.add(a=1.0+1e-2*random.uniform(-1,1), f=0.3*random.uniform(-1,1))
            sim2.add(a=1.0+1e-2*random.uniform(-1,1), f=0.3*random.uniform(-1,1))


        sim.integrate(10)
        sim2.integrate(10)

        self.assertNotEqual(sim.particles[2].x, sim2.particles[2].x)
        self.assertNotEqual(sim.particles[1].x, sim2.particles[1].x)

    def test_tponly_encounter1_mercurius(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-2,a=1)
        sim.add(m=1e-2,a=2)
        sim.move_to_com()
        sim.N_active = sim.N
        sim.N_active = sim.N
        sim.integrator = "mercurius"
        sim.ri_mercurius.tponly_encounter = 1
        sim.dt = 0.1
        sim2 = sim.copy()
        
        random.seed(10)
        for i in range(10):
            sim.add(a=1.0+1e-2*random.uniform(-1,1), f=0.3*random.uniform(-1,1))
            sim2.add(a=1.0+1e-2*random.uniform(-1,1), f=0.3*random.uniform(-1,1))


        sim.integrate(10)
        sim2.integrate(10)

        self.assertEqual(sim.particles[2].x, sim2.particles[2].x)
        self.assertEqual(sim.particles[1].x, sim2.particles[1].x)



if __name__ == "__main__":
    unittest.main()

