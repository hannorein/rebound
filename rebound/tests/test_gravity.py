import rebound
import unittest
import warnings

class TestGravity(unittest.TestCase):
    
    def test_testparticle_0(self):
        sim = rebound.Simulation()
        sim.testparticle_type = 0
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.integrate(10.)
            self.assertEqual(1,len(w))
        x0 = sim.particles[0].x
        self.assertEqual(x0, 0.)

    def test_testparticle_1(self):
        sim = rebound.Simulation()
        sim.testparticle_type = 1
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        sim.integrate(10.)
        x0 = sim.particles[0].x
        self.assertNotEqual(x0, 0.)
    
    def test_testparticle_comp_0(self):
        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.testparticle_type = 0
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.integrate(10.)
            self.assertEqual(1,len(w))
        x0 = sim.particles[0].x
        self.assertEqual(x0, 0.)

    def test_testparticle_comp_1(self):
        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.testparticle_type = 1
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        sim.integrate(10.)
        x0 = sim.particles[0].x
        self.assertNotEqual(x0, 0.)




    def test_testparticle_whfast_comp_0(self):
        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.integrator = "whfast"
        sim.testparticle_type = 0
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.N_active = 1
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.integrate(10.)
            self.assertEqual(1,len(w))
        x0 = sim.particles[0].x
        # TODO
        # Currently fails. WHFAST evolves COM, but should only include star
        #self.assertEqual(x0, 0.)

    def test_testparticle_whfast_comp_1(self):
        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.integrator = "whfast"
        sim.dt = 1e-4
        sim.testparticle_type = 1
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.add(m=1e-3, a=1.4)
        sim.N_active = 2
        sim.integrate(10.)
        x0 = sim.particles[0].x
        x1 = sim.particles[1].x
        self.assertNotEqual(x0, 0.)

        sim = rebound.Simulation()
        sim.gravity = "compensated"
        sim.testparticle_type = 1
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.add(m=1e-3, a=1.4)
        sim.N_active = 2
        sim.integrate(10.)
        x1ias = sim.particles[1].x
        self.assertAlmostEqual(x1ias, x1,delta=1e-9)
    
    def test_tree_duplicate_particle(self):
        sim = rebound.Simulation()
        sim.configure_box(10)
        sim.gravity = "tree"
        sim.add(m=1.)
        with self.assertRaises(RuntimeError):
            # Cannot add two particles on top of each other
            sim.add(m=1.)



if __name__ == "__main__":
    unittest.main()
