import rebound
import unittest
import os

class TestSimulationConstructor(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)

    def tearDown(self):
        self.sim = None

class TestSimulationBinaryCheckpoints(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.t = 1.246
        self.sim.add(m=1.)
        self.sim.add(m=1e-3, a=1., e=0.01, omega=0.02, M=0.04, inc=0.1)
    
    def test_status(self):
        self.sim.status()
    
    def test_N(self):
        self.assertEqual(self.sim.N, 2)
    
    def test_N_real(self):
        self.assertEqual(self.sim.N_real, 2)
    
    def test_integrator(self):
        self.assertEqual(self.sim.integrator, "ias15")
        self.sim.integrator = "whfast"
        self.assertEqual(self.sim.integrator, "whfast")
        self.sim.integrator = 1
        self.assertEqual(self.sim.integrator, "whfast")
    

    def tearDown(self):
        self.sim = None

    def test_nofile(self):
        with self.assertRaises(ValueError):
            sim2 = rebound.Simulation.from_file("doesnotexist.bin")


    def test_status(self):
        self.sim.save("bintest.bin")
        sim2 = rebound.Simulation.from_file("bintest.bin")
        self.assertEqual(self.sim.particles[1].x, sim2.particles[1].x)
        self.assertEqual(self.sim.particles[1].vx, sim2.particles[1].vx)
        self.assertEqual(self.sim.t, sim2.t)
        self.assertEqual(self.sim.N, sim2.N)
        self.assertEqual(self.sim.integrator, sim2.integrator)
        os.remove("bintest.bin")
    
