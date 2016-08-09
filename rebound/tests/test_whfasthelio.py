import rebound
import unittest
import math
import rebound.data
import warnings

class TestIntegratorWHFastHelio(unittest.TestCase):
    def test_whfasthelio_veryhyperbolic(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=0.,x=1.,vy=100000.)
        sim.integrator = "whfasthelio"
        sim.dt = 1.234567
        sim.step()
        y = sim.particles[1].y
        ys = 1.234567*100000.
        self.assertAlmostEqual((y-ys)/ys, 0., delta=1e-10)

    def test_whfasthelio_hyperbolic(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=-1.,e=2.5)
        sim.integrator = "whfasthelio"
        e0 = sim.calculate_energy()
        yr = -sim.particles[1].P
        sim.dt = 0.00512*yr
        sim.integrate(1e2*yr)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-7)
    
    def test_whfasthelio(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.,e=.1)
        sim.add(m=1e-3, a=3.,e=0.1)
        sim.integrator = "whfasthelio"
        jupyr = 2.*math.pi
        sim.dt = 0.005123*jupyr
        e0 = sim.calculate_energy()
        sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-7)
    
    def test_whfasthelio_cor(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.,e=.1)
        sim.add(m=1e-3, a=3.,e=0.1)
        sim.integrator = "whfasthelio"
        sim.ri_whfasthelio.corrector=11
        jupyr = 2.*math.pi
        sim.dt = 0.05123*jupyr
        e0 = sim.calculate_energy()
        sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-7)
    
    def test_whfasthelioi_cor_nosafemode(self):
        sim = rebound.Simulation()
        sim.ri_whfasthelio.safe_mode=0
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.,e=.1)
        sim.add(m=1e-3, a=3.,e=0.1)
        sim.integrator = "whfasthelio"
        sim.ri_whfasthelio.corrector=11
        jupyr = 2.*math.pi
        sim.dt = 0.05123*jupyr
        e0 = sim.calculate_energy()
        sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-7)
    
    
    def test_whfasthelio_orderdoesnotmatter(self):
        jupyr = 2.*math.pi

        sim1 = rebound.Simulation()
        sim1.add(m=1.)
        sim1.add(m=1e-3, a=1.,e=0.1,primary=sim1.particles[0])
        sim1.add(m=1e-3, a=3.,e=0.1,primary=sim1.particles[0])
        sim1.integrator = "whfasthelio"
        sim1.dt = 0.005123*jupyr
        sim1.integrate(1e0*jupyr)
        sim2 = rebound.Simulation()
        sim2.add(m=1.)
        sim2.add(m=1e-3, a=3.,e=0.1,primary=sim2.particles[0])
        sim2.add(m=1e-3, a=1.,e=0.1,primary=sim2.particles[0])
        sim2.integrator = "whfasthelio"
        sim2.dt = 0.005123*jupyr
        sim2.integrate(1e0*jupyr)
        
        self.assertAlmostEqual(sim1.particles[1].x, sim2.particles[2].x, delta=1e-15)
        self.assertAlmostEqual(sim1.particles[1].y, sim2.particles[2].y, delta=1e-15)
        self.assertAlmostEqual(sim1.particles[2].x, sim2.particles[1].x, delta=1e-15)
        self.assertAlmostEqual(sim1.particles[2].y, sim2.particles[1].y, delta=1e-15)
        
    
    def test_whfasthelio_nosafemode(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.,e=.1)
        sim.add(m=1e-3, a=3.,e=0.1)
        sim.integrator = "whfasthelio"
        sim.ri_whfasthelio.safe_mode = 0
        jupyr = 2.*math.pi
        sim.dt = 0.005123*jupyr
        e0 = sim.calculate_energy()
        sim.integrate(1e3*jupyr)
        self.assertNotEqual(e0,0.)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-7)

if __name__ == "__main__":
    unittest.main()
