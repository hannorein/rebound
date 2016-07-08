import rebound
import warnings
import unittest
import os
import math
import sys
from ctypes import c_uint32

class TestSimulation(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.t = 1.246
        self.sim.add(m=1.)
        self.sim.add(m=1e-3, a=1., e=0.01, omega=0.02, M=0.04, inc=0.1)
    
    def tearDown(self):
        self.sim = None

    def test_status(self):
        sys.stdout = open(os.devnull, 'w')
        self.sim.status()
        sys.stdout.close()
        sys.stdout = sys.__stdout__
    
    def test_escape(self):
        self.sim.exit_max_distance = 0.1
        with self.assertRaises(rebound.Escape):
            self.sim.integrate(1.)

    def test_encounter(self):
        self.sim.exit_min_distance = 1.
        with self.assertRaises(rebound.Encounter):
            self.sim.integrate(1.)
    
    def test_removeall(self):
        del self.sim.particles
        self.assertEqual(self.sim.N,0)
    
    def test_remove_too_many(self):
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            self.sim.remove(0)
            self.sim.remove(0)
            self.assertEqual(1,len(w))
        with self.assertRaises(RuntimeError):
            self.sim.remove(0)
    
    def test_remove_variational(self):
        v = self.sim.add_variation()
        with self.assertRaises(RuntimeError):
            self.sim.remove(0)
    
    def test_remove(self):
        self.sim.remove(1)
        self.assertEqual(self.sim.N,1)
    
    def test_remove_keepsorted(self):
        self.sim.remove(1,keepSorted=0)
        self.assertEqual(self.sim.N,1)
    
    def test_removehash(self):
        self.sim.add(m=1e-3, a=1., e=0.01, omega=0.02, M=0.04, inc=0.1)
        self.sim.particles[-1].hash = 99
        self.sim.remove(hash=99)
        self.assertEqual(self.sim.N,2)
        with self.assertRaises(RuntimeError):
            self.sim.remove(hash=99)
        with self.assertRaises(RuntimeError):
            self.sim.remove(hash=99)
        with self.assertRaises(RuntimeError):
            self.sim.remove(hash=-99334)
    
    def test_ascii(self):
        a = self.sim.particles_ascii()
        sim = rebound.Simulation()
        sim.add_particles_ascii(a)
        for i in range(self.sim.N):
            self.assertAlmostEqual(self.sim.particles[i].m,sim.particles[i].m,delta=1e-7)
            self.assertAlmostEqual(self.sim.particles[i].x,sim.particles[i].x,delta=1e-7)
            self.assertAlmostEqual(self.sim.particles[i].vy,sim.particles[i].vy,delta=1e-7)
    
    def test_configure_ghostboxes(self):
        self.sim.configure_ghostboxes(1,1,1)
   
    def test_step(self):
        self.sim.step()
        self.assertNotEqual(self.sim.t, 1.246)

    def test_configure_box(self):
        self.assertEqual(self.sim.root_size,-1.)
        self.sim.configure_box(100.,1,1,1)
        self.assertEqual(self.sim.root_size,100.)
    
    def test_calculate_orbits(self):
        orbits = self.sim.calculate_orbits()
        self.assertAlmostEqual(orbits[0].a,1.,delta=1e-15)
        self.assertAlmostEqual(orbits[0].e,0.01,delta=1e-15)
        self.assertAlmostEqual(orbits[0].omega,0.02,delta=1e-12)
        self.assertAlmostEqual(orbits[0].inc,0.1,delta=1e-15)
        orbits = self.sim.calculate_orbits(heliocentric=True)
        self.assertAlmostEqual(orbits[0].a,1.,delta=1e-15)
        self.assertAlmostEqual(orbits[0].e,0.01,delta=1e-15)
        self.assertAlmostEqual(orbits[0].omega,0.02,delta=1e-12)
        self.assertAlmostEqual(orbits[0].inc,0.1,delta=1e-15)
        orbits = self.sim.calculate_orbits(barycentric=True)
        self.assertAlmostEqual(orbits[0].a,1.,delta=1e-2)
        
    def test_com(self):
        self.sim.move_to_com()
        com = self.sim.calculate_com()
        self.assertAlmostEqual(com.x, 0., delta=1e-15)
    
    def test_com_range(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1., x=2.)
        sim.add(m=2., x=5.)
        com = sim.calculate_com(first=1)
        self.assertAlmostEqual(com.x, 4., delta=1e-15)
        com = sim.calculate_com(last=2)
        self.assertAlmostEqual(com.x, 1., delta=1e-15)
        com = sim.calculate_com(first=1,last=2)
        self.assertAlmostEqual(com.x, 2., delta=1e-15)
        com = sim.calculate_com(first=4, last=-3)
        self.assertAlmostEqual(com.x, 0., delta=1e-15)
    
    def test_jacobi_com(self):
        sim = rebound.Simulation()
        sim.add(m=1., x=1.)
        sim.add(m=1., x=3.)
        sim.add(m=2., x=5.)
        com = sim.particles[1].jacobi_com
        self.assertAlmostEqual(com.x, 1., delta=1e-15)
        com = sim.particles[2].jacobi_com
        self.assertAlmostEqual(com.x, 2., delta=1e-15)
        com = sim.particles[0].jacobi_com
        self.assertAlmostEqual(com.x, 0., delta=1e-15)

    def test_init_megno(self):
        self.sim.init_megno()
        self.assertEqual(self.sim.N,4)
        self.assertEqual(self.sim.N_real,2)
        self.assertEqual(self.sim.calculate_megno(),0)
        self.assertTrue(math.isnan(self.sim.calculate_lyapunov()))
        
    def test_calculate_energy(self):
        self.sim.move_to_com()
        energy = self.sim.calculate_energy()
        self.assertAlmostEqual(energy, -0.5e-3, delta=1e-14)

    def test_calculate_angular_momentum(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-3, a=1., inc=0.3, Omega=0.5)
        sim.add(m=1.e-3, a=3., inc=0.2, Omega = -0.8)
        L0 = sim.calculate_angular_momentum()
        sim.integrate(1.)
        Lf = sim.calculate_angular_momentum()
        for i in range(3):
            self.assertAlmostEqual(abs((Lf[i]-L0[i])/L0[i]), 0., delta=1e-15)

    def test_additional_forces(self):
        def af(sim):
            sim.contents.particles[0].hash = 5
            pass
        self.sim.additional_forces = af
        self.sim.integrate(.1)
        self.assertEqual(self.sim.particles[0].hash.value,5)
        with self.assertRaises(AttributeError):
            self.sim.additional_forces

    def test_post_timestep_modifications(self):
        def ptm(sim):
            sim.contents.particles[0].hash = 6
            pass
        self.sim.post_timestep_modifications = ptm
        self.sim.integrate(.1)
        self.assertEqual(self.sim.particles[0].hash.value,6)
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
    
    def test_checkpoint_ias15_pointers(self):
        self.sim.integrate(1.)
        self.sim.save("bintest.bin")
        self.sim.integrate(5.)
        sim2 = rebound.Simulation.from_file("bintest.bin")
        sim2.integrate(5.)
        self.assertEqual(self.sim.particles[1].x, sim2.particles[1].x)
        self.assertEqual(self.sim.particles[1].vx, sim2.particles[1].vx)
        self.assertEqual(self.sim.t, sim2.t)
        self.assertEqual(self.sim.N, sim2.N)
        self.assertEqual(self.sim.integrator, sim2.integrator)
        os.remove("bintest.bin")
    
    def test_checkpoint_warnings(self):
        with open("bintest2.bin", "a") as f:
            for i in range(1000):
                f.write("gibberish")
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            sim2 = rebound.Simulation.from_file("bintest2.bin")
            self.assertEqual(4,len(w))
        os.remove("bintest2.bin")
    
    
class TestSimulationCollisions(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.gravity = "none"
        self.sim.collision = "direct"
        self.sim.integrator = "leapfrog"
        self.sim.G = 0.0
        self.sim.dt = 0.01
    
    def tearDown(self):
        self.sim = None

    def test_coefficient_of_restitution(self):
        self.sim.add(m=1.,x=-1,vx=1.,r=0.5)
        self.sim.add(m=1.,x=1,vx=-1.,r=0.5)
        energy_initial = self.sim.calculate_energy()
        def coef(sim,vrel):
            return 0.5
        self.sim.coefficient_of_restitution = coef
        self.sim.integrate(1.)
        energy_final = self.sim.calculate_energy()
        self.assertAlmostEqual(energy_final, 0.25*energy_initial,delta=1e-15)

    def test_direct(self):
        self.sim.add(m=1.,x=-1,vx=1.,r=0.5)
        self.sim.add(m=1.,x=1,vx=-1.,r=0.5)
        self.sim.integrate(1.)
        self.assertAlmostEqual(self.sim.particles[0].x,-1,delta=1e-15)

    def test_tree(self):
        self.sim.configure_box(10)
        self.sim.collision = "tree"
        self.sim.add(m=1.,x=-1,vx=1.,r=0.5)
        self.sim.add(m=1.,x=1,vx=-1.,r=0.5)
        self.sim.integrate(1.)
        self.assertAlmostEqual(self.sim.particles[0].x,-1,delta=1e-15)

    
    
if __name__ == "__main__":
    unittest.main()
