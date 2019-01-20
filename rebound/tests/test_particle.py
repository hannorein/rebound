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

    def test_masses(self):
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3, a=1.)
        self.sim.add(m=1.e-6, a=2.)
        self.sim.add(m=0., a=3.)
        ps = self.sim.particles
        self.assertAlmostEqual(ps[0].m,1.,delta=1e-15)
        self.assertAlmostEqual(ps[1].m,1.e-3,delta=1e-15)
        self.assertAlmostEqual(ps[2].m,1.e-6,delta=1e-15)
        self.assertAlmostEqual(ps[3].m,0.,delta=1e-15)

    def test_jacobi_masses(self):
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3, a=1., jacobi_masses=True)
        self.sim.add(m=1.e-6, a=2., jacobi_masses=True)
        self.sim.add(m=0., a=3., jacobi_masses=True)
        ps = self.sim.particles
        self.assertAlmostEqual(ps[0].m,1.,delta=1e-15)
        self.assertAlmostEqual(ps[1].m,1.e-3,delta=1e-15)
        self.assertAlmostEqual(ps[2].m,1.e-6,delta=1e-15)
        self.assertAlmostEqual(ps[3].m,0.,delta=1e-15)
        self.assertAlmostEqual(ps[1].vy,1.000499875062461,delta=1e-15)
        self.assertAlmostEqual(ps[2].vy,0.7081066347613374,delta=1e-15)
        self.assertAlmostEqual(ps[3].vy,0.5783504759643414,delta=1e-15)
        o = self.sim.calculate_orbits(jacobi_masses=True)
        self.assertAlmostEqual(o[0].a,1.,delta=1e-15)
        self.assertAlmostEqual(o[1].a,2.,delta=1e-15)
        self.assertAlmostEqual(o[2].a,3.,delta=1e-15)
        self.assertAlmostEqual(o[0].e,0.,delta=1e-15)
        self.assertAlmostEqual(o[1].e,0.,delta=1e-15)
        self.assertAlmostEqual(o[2].e,0.,delta=1e-15)

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
   
    def test_sim_calculate_orbits(self):
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3, a=1.,e=0.2,inc=0.3)
        self.sim.add(m=1.e-3, a=2.,e=0.2,inc=0.3)
        ps = self.sim.particles
        self.assertAlmostEqual(ps[1].a,1.,delta=1e-15)
        self.assertAlmostEqual(ps[1].e,0.2,delta=1e-15)
        self.assertAlmostEqual(ps[1].inc,0.3,delta=1e-15)
        self.assertAlmostEqual(ps[2].a,2.,delta=1e-15)
        self.assertAlmostEqual(ps[2].e,0.2,delta=1e-15)
        self.assertAlmostEqual(ps[2].inc,0.3,delta=1e-15)

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
        self.assertEqual(p3.x,p1.x-p2.x)
        self.assertEqual(p3.y,p1.y-p2.y)
        self.assertEqual(p3.z,p1.z-p2.z)
        self.assertEqual(p3.vx,p1.vx-p2.vx)
        self.assertEqual(p3.vy,p1.vy-p2.vy)
        self.assertEqual(p3.vz,p1.vz-p2.vz)
        self.assertEqual(p3.m,p1.m-p2.m)

    def test_add(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = rebound.Particle(m=1.4,x=1.6,vy=1.8)
        p3 = p1 + p2
        self.assertEqual(p3.x,p1.x+p2.x)
        self.assertEqual(p3.y,p1.y+p2.y)
        self.assertEqual(p3.z,p1.z+p2.z)
        self.assertEqual(p3.vx,p1.vx+p2.vx)
        self.assertEqual(p3.vy,p1.vy+p2.vy)
        self.assertEqual(p3.vz,p1.vz+p2.vz)
        self.assertEqual(p3.m,p1.m+p2.m)

    def test_mul(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = 2.*p1
        self.assertEqual(p2.x,2.*p1.x)
        self.assertEqual(p2.y,2.*p1.y)
        self.assertEqual(p2.z,2.*p1.z)
        self.assertEqual(p2.vx,2.*p1.vx)
        self.assertEqual(p2.vy,2.*p1.vy)
        self.assertEqual(p2.vz,2.*p1.vz)
        self.assertEqual(p2.m,2.*p1.m)
    
    def test_div(self):
        p1 = rebound.Particle(m=1.2,x=1.4,vy=1.8)
        p2 = p1/2.
        self.assertEqual(p2.x,p1.x/2.)
        self.assertEqual(p2.y,p1.y/2.)
        self.assertEqual(p2.z,p1.z/2.)
        self.assertEqual(p2.vx,p1.vx/2.)
        self.assertEqual(p2.vy,p1.vy/2.)
        self.assertEqual(p2.vz,p1.vz/2.)
        self.assertEqual(p2.m,p1.m/2.)
    
    def test_truediv(self):
        p1 = rebound.Particle(m=1.2,x=1.4,vy=1.8)
        p2 = p1/2
        self.assertEqual(p2.x,p1.x/2.)
        self.assertEqual(p2.y,p1.y/2.)
        self.assertEqual(p2.z,p1.z/2.)
        self.assertEqual(p2.vx,p1.vx/2.)
        self.assertEqual(p2.vy,p1.vy/2.)
        self.assertEqual(p2.vz,p1.vz/2.)
        self.assertEqual(p2.m,p1.m/2.)
    
    def test_isub(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p3 = p1.copy()
        p2 = rebound.Particle(m=1.4,x=1.6,vy=1.8)
        p3 -= p2
        self.assertEqual(p3.x,p1.x-p2.x)
        self.assertEqual(p3.y,p1.y-p2.y)
        self.assertEqual(p3.z,p1.z-p2.z)
        self.assertEqual(p3.vx,p1.vx-p2.vx)
        self.assertEqual(p3.vy,p1.vy-p2.vy)
        self.assertEqual(p3.vz,p1.vz-p2.vz)
        self.assertEqual(p3.m,p1.m-p2.m)

    def test_iadd(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p3 = p1.copy()
        p2 = rebound.Particle(m=1.4,x=1.6,vy=1.8)
        p3 += p2
        self.assertEqual(p3.x,p1.x+p2.x)
        self.assertEqual(p3.y,p1.y+p2.y)
        self.assertEqual(p3.z,p1.z+p2.z)
        self.assertEqual(p3.vx,p1.vx+p2.vx)
        self.assertEqual(p3.vy,p1.vy+p2.vy)
        self.assertEqual(p3.vz,p1.vz+p2.vz)
        self.assertEqual(p3.m,p1.m+p2.m)

    def test_imul(self):
        p1 = rebound.Particle(m=1.1,x=1.2,vy=1.3)
        p2 = p1.copy()
        p2 *= 2.
        self.assertEqual(p2.x,2.*p1.x)
        self.assertEqual(p2.y,2.*p1.y)
        self.assertEqual(p2.z,2.*p1.z)
        self.assertEqual(p2.vx,2.*p1.vx)
        self.assertEqual(p2.vy,2.*p1.vy)
        self.assertEqual(p2.vz,2.*p1.vz)
        self.assertEqual(p2.m,2.*p1.m)
    
    def test_idiv(self):
        p1 = rebound.Particle(m=1.2,x=1.4,vy=1.8)
        p2 = p1.copy()
        p2 /=2.
        self.assertEqual(p2.x,p1.x/2.)
        self.assertEqual(p2.y,p1.y/2.)
        self.assertEqual(p2.z,p1.z/2.)
        self.assertEqual(p2.vx,p1.vx/2.)
        self.assertEqual(p2.vy,p1.vy/2.)
        self.assertEqual(p2.vz,p1.vz/2.)
        self.assertEqual(p2.m,p1.m/2.)
    
    def test_itruediv(self):
        p1 = rebound.Particle(m=1.2,x=1.4,vy=1.8)
        p2 = p1.copy()
        p2 /=2
        self.assertEqual(p2.x,p1.x/2.)
        self.assertEqual(p2.y,p1.y/2.)
        self.assertEqual(p2.z,p1.z/2.)
        self.assertEqual(p2.vx,p1.vx/2.)
        self.assertEqual(p2.vy,p1.vy/2.)
        self.assertEqual(p2.vz,p1.vz/2.)
        self.assertEqual(p2.m,p1.m/2.)
    


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
