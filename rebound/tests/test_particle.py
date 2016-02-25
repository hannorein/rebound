import rebound
import unittest
import math 

class TestParticleInSimulation(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_adding(self):
        self.sim.add()
        self.sim.add(m=1.)
        self.sim.add(m=1.,x=1, vy=1)
        self.sim.add(x=2, id=9)
        self.sim.add(x=3, r=1.)
        with self.assertRaises(ValueError):
            self.sim.add(x=4,a=1)

    def test_adding_orbits(self):
        self.sim.add(m=1.)
        self.sim.add(a=1.)
        with self.assertRaises(ValueError):
            self.sim.add(e=0.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=1.,e=0.1,omega=0.1,pomega=0.1)
        self.sim.add(a=2.,e=0.1,inc=0.1,pomega=0.1)
        self.sim.add(a=2.,e=0.1,inc=-0.1,pomega=0.1)
        self.sim.add(a=2.,e=0.1,inc=0.1,theta=0.1)
        self.sim.add(a=2.,e=0.1,inc=-0.1,theta=0.1)
        self.sim.add(a=2.,e=0.1,inc=0.1,l=0.1)
        self.sim.add(a=2.,e=0.1,inc=-0.1,l=0.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=2.,e=0.1,f=0.1,M=0.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=2.,e=0.1,f=0.1,l=0.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=2.,e=0.1,f=0.1,theta=0.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=3.,e=1.)
        with self.assertRaises(ValueError):
            self.sim.add(a=3.,e=1.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=3.,e=-0.1)
        self.sim.add(a=-3.,e=1.4)
        with self.assertRaises(ValueError):
            self.sim.add(a=-3.,e=0.9)
        self.sim.add(a=-3.,e=1.4,f=0.1)
        with self.assertRaises(ValueError):
            self.sim.add(a=-3.,e=1.4,f=3.1)
    
    def test_calculate_orbits(self):
        self.sim.add(m=1.)
        self.sim.add(a=1.)
        p = self.sim.particles
        with self.assertRaises(ValueError):
            p[0].calculate_orbit()
        o = p[1].calculate_orbit()
        self.assertAlmostEqual(o.a,1.,delta=1e-15)
        self.assertAlmostEqual(o.e,0.,delta=1e-15)
        self.assertAlmostEqual(o.f,0.,delta=1e-15)
        self.assertAlmostEqual(o.inc,0.,delta=1e-15)
        
        self.assertAlmostEqual(p[1].a,1.,delta=1e-15)
        self.assertAlmostEqual(p[1].e,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].f,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].inc,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].d,1.,delta=1e-15)
        self.assertAlmostEqual(p[1].v,1.,delta=1e-15)
        self.assertAlmostEqual(p[1].h,1.,delta=1e-15)
        self.assertAlmostEqual(p[1].P,math.pi*2.,delta=1e-15)
        self.assertAlmostEqual(p[1].n,1.,delta=1e-15)
        self.assertAlmostEqual(p[1].omega,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].pomega,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].Omega,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].M,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].l,0.,delta=1e-15)
        self.assertAlmostEqual(p[1].theta,0.,delta=1e-15)
    
    def test_calculate_orbits_errors(self):
        self.sim.add()
        self.sim.add(x=1)
        with self.assertRaises(ValueError):
            self.sim.particles[1].orbit

    def test_calculate_orbits_errors2(self):
        self.sim.add(m=1)
        p1 = rebound.Particle(simulation=self.sim, a=1,m=0.1)
        self.sim.add(p1)
        with self.assertRaises(ValueError):
            self.sim.particles[1].calculate_orbit(primary=p1)



class TestParticleNotInSimulation(unittest.TestCase):
    def test_create(self):
        p1 = rebound.Particle()
        p2 = rebound.Particle(x=1)
        with self.assertRaises(ValueError):
            p3 = rebound.Particle(a=1)
        with self.assertRaises(ValueError):
            p4 = rebound.Particle(p1)

