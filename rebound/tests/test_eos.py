import rebound
import unittest
import math 

class TestEOSn(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1)
        self.sim.add(m=1e-3,a=1,e=0.05,f=0.)
        self.sim.add(m=1e-3,a=1.6,e=0.05,f=1.)
        self.sim.move_to_com()
    
    def tearDown(self):
        self.sim = None

    def _run(self):
        tmax = 100.
        Emax = 0
        E0 = self.sim.calculate_energy()
        while self.sim.t<tmax:
            self.sim.integrate(self.sim.t+1.23456, exact_finish_time=0)
            E1 = self.sim.calculate_energy()
            Emax = max([Emax,abs((E0-E1)/E0)])
        return Emax
    
    def test_lf(self):
        self.sim.dt = 0.01*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "lf"
        self.sim.ri_eos.phi1 = "lf"
        self.sim.ri_eos.n = 16
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=5.1e-7)
    
    def test_lf4(self):
        self.sim.dt = 0.001*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "lf4"
        self.sim.ri_eos.phi1 = "lf4"
        self.sim.ri_eos.n = 2
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=2.6e-11)
    
    def test_lf4_2(self):
        self.sim.dt = 0.001*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "lf4_2"
        self.sim.ri_eos.phi1 = "lf4"
        self.sim.ri_eos.n = 2
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=7.0e-12)
    
    def test_pmlf4(self):
        self.sim.dt = 0.001*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "pmlf4"
        self.sim.ri_eos.phi1 = "lf4"
        self.sim.ri_eos.n = 2
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=7.0e-12)
    
    def test_lf6(self):
        self.sim.dt = 0.01*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "lf6"
        self.sim.ri_eos.phi1 = "lf6"
        self.sim.ri_eos.n = 2
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=8e-14)
    
    def test_lf8(self):
        self.sim.dt = 0.03*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "lf8"
        self.sim.ri_eos.phi1 = "lf8"
        self.sim.ri_eos.n = 1
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=2e-13)
    
    def test_lf8_6_4(self):
        self.sim.dt = 0.03*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "lf8_6_4"
        self.sim.ri_eos.phi1 = "lf8"
        self.sim.ri_eos.n = 1
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=1.4e-13)
    
    def test_plf7_6_4(self):
        self.sim.dt = 0.01*(2.*math.pi)
        self.sim.integrator = "eos"
        self.sim.ri_eos.phi0 = "plf7_6_4"
        self.sim.ri_eos.phi1 = "lf8"
        self.sim.ri_eos.n = 1
        self.sim.ri_eos.safe_mode = 0
        Emax = self._run()
        self.assertAlmostEqual(Emax,0.,delta=1.5e-13)
    
    
if __name__ == "__main__":
    unittest.main()
