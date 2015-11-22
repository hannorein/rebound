import rebound
import unittest
import math

class TestIntegrator(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(self.sim)
        self.sim.move_to_com()
    
    def tearDown(self):
        self.sim = None
    
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
        self.sim.integrator_whfast_safe_mode = 0
        self.sim.integrator_whfast_corrector = 11
        jupyr = 11.86*2.*math.pi
        self.sim.dt = 0.0123*jupyr
        e0 = self.sim.calculate_energy()
        self.sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = self.sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-9)
