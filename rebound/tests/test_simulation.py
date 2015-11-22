import rebound
import unittest
import os

class TestSimulation(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.t = 1.246
        self.sim.add(m=1.)
        self.sim.add(m=1e-3, a=1., e=0.01, omega=0.02, M=0.04, inc=0.1)
    
    def tearDown(self):
        self.sim = None

    def test_status(self):
        self.sim.status()
    
    def test_calculate_orbits(self):
        orbits = self.sim.calculate_orbits()
        print orbits[0]
        
    def test_com(self):
        self.sim.move_to_com()
        com = self.sim.calculate_com()
        self.assertEqual(com.x, 0.)
        
        
    def test_calculate_energy(self):
        self.sim.move_to_com()
        energy = self.sim.calculate_energy()
        self.assertAlmostEqual(energy, -0.5e-3, delta=1e-14)


    def test_additional_forces(self):
        def af(sim):
            sim.contents.particles[0].id = 5
            pass
        self.sim.additional_forces = af
        self.sim.integrate(.1)
        self.assertEqual(self.sim.particles[0].id,5)
        with self.assertRaises(AttributeError):
            self.sim.additional_forces

    def test_post_timestep_modifications(self):
        def ptm(sim):
            sim.contents.particles[0].id = 6
            pass
        self.sim.post_timestep_modifications = ptm
        self.sim.integrate(.1)
        self.assertEqual(self.sim.particles[0].id,6)
        with self.assertRaises(AttributeError):
            self.sim.post_timestep_modifications


    
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
        with self.assertRaises(ValueError):
            self.sim.integrator = "bogusintegrator"

    def test_boundaries(self):
        self.sim.boundary = "open"
        self.assertEqual(self.sim.boundary, "open")
        self.sim.boundary = 8
        self.assertEqual(self.sim.boundary, 8)
        with self.assertRaises(ValueError):
            self.sim.boundary = "bogusboundary"

    
    def test_gravity(self):
        self.sim.gravity = "tree"
        self.assertEqual(self.sim.gravity, "tree")
        self.sim.gravity = 8
        self.assertEqual(self.sim.gravity, 8)
        with self.assertRaises(ValueError):
            self.sim.gravity = "bogusgravity"
    
    def test_collision(self):
        self.sim.collision = "tree"
        self.assertEqual(self.sim.collision, "tree")
        self.sim.collision = 8
        self.assertEqual(self.sim.collision, 8)
        with self.assertRaises(ValueError):
            self.sim.collision = "boguscollision"


    def test_nofile(self):
        with self.assertRaises(ValueError):
            sim2 = rebound.Simulation.from_file("doesnotexist.bin")


    def test_checkpoint(self):
        self.sim.save("bintest.bin")
        sim2 = rebound.Simulation.from_file("bintest.bin")
        self.assertEqual(self.sim.particles[1].x, sim2.particles[1].x)
        self.assertEqual(self.sim.particles[1].vx, sim2.particles[1].vx)
        self.assertEqual(self.sim.t, sim2.t)
        self.assertEqual(self.sim.N, sim2.N)
        self.assertEqual(self.sim.integrator, sim2.integrator)
        os.remove("bintest.bin")
    
