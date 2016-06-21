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
        self.sim.add(x=2)
        self.sim.add(x=3, r=1.)
        with self.assertRaises(ValueError):
            self.sim.add(x=4,a=1)
        p4 = self.sim.particles[4]
        ind = p4.index
        self.assertEqual(ind,4)

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
        self.sim.add(P=2.,e=0.1,inc=-2.1,pomega=0.1)
        self.sim.add(P=2.,e=0.1,inc=-2.1,pomega=0.1,f=0.2)
        self.sim.add(P=2.,e=0.1,inc=-2.1,pomega=0.1,T=0.2)
        self.sim.add(P=2.,e=0.1,inc=-2.1,pomega=0.1,theta=0.2)
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
            self.sim.add(a=3.,P=1.1)
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
        string = o.__str__()
        self.assertGreater(len(string),20)
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
        self.assertAlmostEqual(p[1].T,0.,delta=1e-15)
    
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
    
    def test_calculate_orbits_errors3(self):
        p1 = rebound.Particle(m=1.,x=1.,vy=0.4)
        p2 = rebound.Particle(m=1.,x=4.,vy=2.4)
        with self.assertRaises(ValueError):
            p2.calculate_orbit()
        with self.assertRaises(ValueError):
            p2.calculate_orbit(primary=p1)
        p2.calculate_orbit(primary=p1,G=1.)
        

class TestParticleOperators(unittest.TestCase):
    def test_sub(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = rebound.Particle(m=1.4,x=1.6,vy=1.8)
        p3 = p1 - p2
        self.assertAlmostEqual(p3.m,-0.3,delta=1e-15)
        self.assertAlmostEqual(p3.x,-0.4,delta=1e-15)
        self.assertAlmostEqual(p3.vy,-0.5,delta=1e-15)

    def test_add(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = rebound.Particle(m=1.4,x=1.6,vy=1.8)
        p3 = p1 + p2
        self.assertAlmostEqual(p3.m,2.5,delta=1e-15)
        self.assertAlmostEqual(p3.x,2.8,delta=1e-15)
        self.assertAlmostEqual(p3.vy,3.1,delta=1e-15)

    def test_mul(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = 2.*p1
        self.assertAlmostEqual(p2.m,2.2,delta=1e-15)
        self.assertAlmostEqual(p2.x,2.4,delta=1e-15)
        self.assertAlmostEqual(p2.vy,2.6,delta=1e-15)
    
    def test_irmul(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = p1*2.
        self.assertAlmostEqual(p2.m,2.2,delta=1e-15)
        self.assertAlmostEqual(p2.x,2.4,delta=1e-15)
        self.assertAlmostEqual(p2.vy,2.6,delta=1e-15)
    
    def test_div(self):
        p1 = rebound.Particle(m=1.2,x=1.4,vy=1.8)
        p2 = p1/2.
        self.assertAlmostEqual(p2.m,0.6,delta=1e-15)
        self.assertAlmostEqual(p2.x,0.7,delta=1e-15)
        self.assertAlmostEqual(p2.vy,0.9,delta=1e-15)
    


class TestParticleCopy(unittest.TestCase):
    def test_copy(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.,e=0.2,a=1)
        pc = sim.particles[1].copy()
        self.assertEqual(sim.particles[1].x,pc.x)
        self.assertEqual(sim.particles[1].vx,pc.vx)
        self.assertEqual(sim.particles[1].m,pc.m)
        sim.particles[1].m=0.01
        self.assertNotEqual(sim.particles[1].m,pc.m)
    def test_copy2(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        p1 = rebound.Particle(simulation=sim,m=1.,e=0.2,a=1)
        p2 = rebound.Particle(p1)
        p1.m=2.
        p2.m=3.
        sim.add(p1)
        sim.add(p2)
        self.assertEqual(p1.m,2.)
        self.assertEqual(p2.m,3.)
        self.assertEqual(sim.particles[1].m,2.)
        self.assertEqual(sim.particles[2].m,3.)
        

class TestParticleNotInSimulation(unittest.TestCase):
    def test_create(self):
        p1 = rebound.Particle()
        p2 = rebound.Particle(x=1)
        with self.assertRaises(ValueError):
            p3 = rebound.Particle(a=1)

if __name__ == "__main__":
    unittest.main()
