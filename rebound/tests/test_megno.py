import rebound
import unittest
import math 

class TestMegno(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_ias15(self):
        self.sim.integrator = "ias15"
        self.sim.add(m=1)
        self.sim.add(m=1e-3,a=1.5,e=0.1,inc=0.1)
        self.sim.init_megno()
        self.sim.integrate(1000)
        self.assertAlmostEqual(self.sim.calculate_megno(),2.,delta=2e-1)
        self.assertAlmostEqual(self.sim.calculate_lyapunov(),0.,delta=1e-3)

    def test_whfast(self):
        self.sim.integrator = "whfast"
        self.sim.add(m=1)
        self.sim.add(m=1e-3,a=1.5,e=0.1,inc=0.1)
        self.sim.init_megno()
        self.sim.integrate(1000)
        self.assertAlmostEqual(self.sim.calculate_megno(),2.,delta=2e-1)
        self.assertAlmostEqual(self.sim.calculate_lyapunov(),0.,delta=1e-3)


    
if __name__ == "__main__":
    unittest.main()
