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
        self.sim.add(m=1.e-3, a=15., e=0.1, inc=0.1)
        self.sim.init_megno(seed=0)
        self.sim.integrate(10000.)
        self.assertAlmostEqual(self.sim.calculate_megno(),2.,delta=2e-1)
        self.assertAlmostEqual(self.sim.calculate_lyapunov(),0.,delta=1e-3)

    def test_whfast(self):
        self.sim.integrator = "whfast"
        self.sim.add(m=1)
        self.sim.add(m=1e-3,a=1.5,e=0.1,inc=0.1)
        self.sim.add(m=1.e-3, a=15., e=0.1, inc=0.1)
        self.sim.init_megno(seed=0)
        self.sim.dt = self.sim.particles[1].P*0.07
        self.sim.integrate(10000.)
        self.assertAlmostEqual(self.sim.calculate_megno(),2.,delta=2e-1)
        self.assertAlmostEqual(self.sim.calculate_lyapunov(),0.,delta=1e-3)

    def test_whfast_close_regular(self):
        self.sim = rebound.Simulation()
        self.sim.integrator = "whfast"
        self.sim.G = 4*math.pi**2
        self.sim.add(m=1.0, x=-4.7298443749900997e-07, y=2.3452237398515157e-06, z=-2.8779986923394045e-08, vx=-9.994358883493113e-06, vy=-2.8416720717989476e-06, vz=1.5927607755978474e-07)
        self.sim.add(m=4.544728982917554e-07, x=0.9818659845399763, y=0.19089291910956852, z=-0.011300485084842085, vx=-1.1951384526938573, vy=6.164295939653557, vz=0.16345519841419276)
        self.sim.add(m=1.7367161546229798e-06, x=-0.07616371745657412, y=-1.2648949778387153, z=0.01873838700271951, vx=5.516651885197617, vy=-0.34889260830083374, vz=-0.13577609379397024)
        self.sim.add(m=2.1454312223049496e-07, x=0.741238938912468, y=-1.0963570106709737, z=0.006397276680495443, vx=4.459049824337919, vy=3.011488098634852, vz=0.010452444973871466)
        self.sim.init_megno(seed=0)
        self.sim.dt = 0.034641008279678746
        self.sim.integrate(10000.)
        self.assertAlmostEqual(self.sim.calculate_megno(),2.,delta=2e-1)
        self.assertAlmostEqual(self.sim.calculate_lyapunov(),0.,delta=1e-3)

    def test_chaotic(self):
        self.sim = rebound.Simulation()
        self.sim.integrator = "ias15"
        self.sim.add(m=1.)
        self.sim.add(m=1.e-4, P=1.)
        self.sim.add(m=1.e-4, P=1.17)
        self.sim.init_megno(seed=0)
        self.sim.move_to_com()
        self.sim.integrate(1000.)
        self.megnoIAS = self.sim.calculate_megno()
        self.sim = rebound.Simulation()
        self.sim.integrator = "whfast"
        self.sim.add(m=1.)
        self.sim.add(m=1.e-4, P=1.)
        self.sim.add(m=1.e-4, P=1.17)
        self.sim.init_megno(seed=0)
        self.sim.move_to_com()
        self.sim.integrate(1000)
        self.megnoWHFast = self.sim.calculate_megno()
        self.assertAlmostEqual(abs((self.megnoIAS-self.megnoWHFast)/self.megnoIAS), 0., delta=0.3)

if __name__ == "__main__":
    unittest.main()
