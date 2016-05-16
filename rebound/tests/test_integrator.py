import rebound
import unittest
import math
import rebound.data

class TestIntegrator2(unittest.TestCase):
    def test_whfast_verylargedt(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.move_to_com()
        sim.integrator = "whfast"
        yr = sim.particles[1].P
        sim.dt = 4.56*yr
        x0 = sim.particles[1].x
        sim.integrate(1e3*yr)
        x1 = sim.particles[1].x
        self.assertAlmostEqual(x0, x1, delta=1e-12)
    
    def test_wh_verylargedt(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.integrator = "wh"
        yr = sim.particles[1].P
        sim.dt = 4.56*yr
        x0 = sim.particles[1].x
        sim.integrate(1e3*yr)
        x1 = sim.particles[1].x
        self.assertAlmostEqual(x0, x1, delta=1e-12)
    
    def test_whfast_hyperbolic(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=-1.,e=2.5)
        sim.integrator = "whfast"
        x0 = sim.calculate_energy()
        yr = -sim.particles[1].P
        sim.dt = 0.12*yr
        sim.integrate(1e2*yr)
        x1 = sim.calculate_energy()
        self.assertAlmostEqual(x0, x1, delta=1e-14)
    
    def test_whfast_verylargedt_hyperbolic(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=-1.,e=2.5)
        sim.integrator = "whfast"
        x0 = sim.calculate_energy()
        yr = -sim.particles[1].P
        sim.dt = 4.56*yr
        sim.integrate(1e3*yr)
        x1 = sim.calculate_energy()
        self.assertAlmostEqual(x0, x1, delta=1e-14)


class TestIntegrator(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(self.sim)
        self.sim.move_to_com()
    
    def tearDown(self):
        self.sim = None
    
    def test_ias15_globaloff(self):
        self.sim.integrator = "ias15"
        self.sim.ri_ias15.epsilon_global = 0
        jupyr = 11.86*2.*math.pi
        e0 = self.sim.calculate_energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3*jupyr)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-14)

    def test_ias15_small_initial_dt(self):
        self.sim.integrator = "ias15"
        jupyr = 11.86*2.*math.pi
        self.sim.dt = jupyr*1e-7
        e0 = self.sim.calculate_energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(self.sim.dt*1.001)
        self.sim.dt = jupyr
        self.sim.integrate(1e3*jupyr)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-14)

    def test_ias15(self):
        self.sim.integrator = "ias15"
        jupyr = 11.86*2.*math.pi
        e0 = self.sim.calculate_energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3*jupyr)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-14)
        # Longer tests:
        #self.sim.integrate(1e4*jupyr)
        #e1 = self.sim.calculate_energy()
        #self.assertLess(math.fabs((e0-e1)/e1),10**13.5)
    
    def test_ias15_compensated(self):
        self.sim.integrator = "ias15"
        self.sim.gravity = "compensated"
        jupyr = 11.86*2.*math.pi
        e0 = self.sim.calculate_energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3*jupyr)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-14)
    
    def test_wh(self):
        self.sim.integrator = "wh"
        self.sim.move_to_com()
        e0 = self.sim.calculate_energy()
        # Move to heliocentric frame
        sun = self.sim.particles[0].copy()
        for p in self.sim.particles:
            m = p.m
            p -= sun
            p.m = m
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.123*jupyr
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        self.sim.move_to_com()
        e1 = self.sim.calculate_energy()
        # Something wrong here with the energy conservations! TODO!
        self.assertLess(math.fabs((e0-e1)/e1),1e-2)
    
    def test_whfast_largedt(self):
        self.sim.integrator = "whfast"
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.123*jupyr
        e0 = self.sim.calculate_energy()
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-4)
    
    def test_whfast_smalldt(self):
        self.sim.integrator = "whfast"
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.0123*jupyr
        e0 = self.sim.calculate_energy()
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-6)
    
    def test_whfast_smalldt_compensated(self):
        self.sim.integrator = "whfast"
        self.sim.gravity = "compensated"
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.0123*jupyr
        e0 = self.sim.calculate_energy()
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-6)
    
    def test_whfast_verysmalldt(self):
        self.sim.integrator = "whfast"
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.00123*jupyr
        e0 = self.sim.calculate_energy()
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-8)
    
    def test_whfast_nosafemode(self):
        self.sim.integrator = "whfast"
        self.sim.ri_whfast.safe_mode = 0
        self.sim.ri_whfast.corrector = 11
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.0123*jupyr
        e0 = self.sim.calculate_energy()
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-9)

if __name__ == "__main__":
    unittest.main()
