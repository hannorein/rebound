import rebound
import math
import unittest
import ctypes

class TestOrbitalElements(unittest.TestCase):
    def test_add_errors(self):
        sim = rebound.Simulation()
        with self.assertRaises(ValueError):
            sim.add(a=1)
        sim.add(m=1.)
        with self.assertRaises(ValueError):
            sim.add(a=1, e=1)
        with self.assertRaises(ValueError):
            sim.add(a=1, e=-2)
        with self.assertRaises(ValueError):
            sim.add(a=1, e=1.2)
        with self.assertRaises(ValueError):
            sim.add(a=-1, e=.2)
        with self.assertRaises(ValueError):
            sim.add(a=-1, e=2.2, f=3.)
    
    def test_orbit_errors(self):
        sim = rebound.Simulation()
        sim.add()
        sim.add(x=1.)
        with self.assertRaises(ValueError):
            a = sim.particles[1].a
        
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-3)
        with self.assertRaises(ValueError):
            a = sim.particles[1].a

    def test_inclined_eccentric(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,l=0.632)
        sim.add(m=1.e-4,a=2.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,f=0.632)
        sim.add(m=1.e-4,a=3.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,theta=0.632)
        sim.add(m=1.e-4,a=4.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,M=0.632)
        sim.add(m=1.e-4,a=5.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,T=0.632)
        sim.add(m=1.e-4,a=6.24,e=0.123,inc=0.14,pomega=0.12,Omega=0.64,l=0.632)
        sim.add(m=1.e-4,a=7.24,e=0.123,inc=0.14,pomega=0.12,Omega=0.64,f=0.632)
        sim.add(m=1.e-4,a=8.24,e=0.123,inc=0.14,pomega=0.12,Omega=0.64,theta=0.632)
        sim.add(m=1.e-4,a=9.24,e=0.123,inc=0.14,pomega=0.12,Omega=0.64,M=0.632)
        sim.add(m=1.e-4,a=10.24,e=0.123,inc=0.14,pomega=0.12,Omega=0.64,T=0.632)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0.123, delta=d)
            self.assertAlmostEqual(p.inc, 0.14, delta=d)
            self.assertAlmostEqual(p.Omega, 0.64, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[1].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[1].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[1].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[2].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[2].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[2].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 3.24, delta=d)
        self.assertAlmostEqual(ps[3].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[3].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[4].a, 4.24, delta=d)
        self.assertAlmostEqual(ps[4].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[4].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[5].a, 5.24, delta=d)
        self.assertAlmostEqual(ps[5].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[5].T, 0.632, delta=d)
        self.assertAlmostEqual(ps[6].a, 6.24, delta=d)
        self.assertAlmostEqual(ps[6].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[6].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[7].a, 7.24, delta=d)
        self.assertAlmostEqual(ps[7].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[7].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[8].a, 8.24, delta=d)
        self.assertAlmostEqual(ps[8].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[8].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[9].a, 9.24, delta=d)
        self.assertAlmostEqual(ps[9].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[9].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[10].a, 10.24, delta=d)
        self.assertAlmostEqual(ps[10].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[10].T, 0.632, delta=d)
    
    def test_planar_eccentric(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.00,e=0.123) # dummy, make sure jacobi com is calculated right for higher indices
        sim.add(m=1.e-4,a=1.24,e=0.123,pomega=0.12,l=0.632)
        sim.add(m=1.e-4,a=2.24,e=0.123,pomega=0.12,f=0.632)
        sim.add(m=1.e-4,a=3.24,e=0.123,pomega=0.12,theta=0.632)
        sim.add(m=1.e-4,a=4.24,e=0.123,pomega=0.12,M=0.632)
        sim.add(m=1.e-4,a=5.24,e=0.123,pomega=0.12,T=0.632)
        sim.add(m=1.e-4,a=6.24,e=0.123,Omega=1.25,pomega=0.12,l=0.632)
        sim.add(m=1.e-4,a=7.24,e=0.123,Omega=1.25,pomega=0.12,f=0.632)
        sim.add(m=1.e-4,a=8.24,e=0.123,Omega=1.25,pomega=0.12,theta=0.632)
        sim.add(m=1.e-4,a=9.24,e=0.123,Omega=1.25,pomega=0.12,M=0.632)
        sim.add(m=1.e-4,a=10.24,e=0.123,Omega=1.25,pomega=0.12,T=6.32)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0.123, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[2].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[2].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[2].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[3].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[3].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[4].a, 3.24, delta=d)
        self.assertAlmostEqual(ps[4].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[4].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[5].a, 4.24, delta=d)
        self.assertAlmostEqual(ps[5].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[5].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[6].a, 5.24, delta=d)
        self.assertAlmostEqual(ps[6].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[6].T, 0.632, delta=d)
        self.assertAlmostEqual(ps[7].a, 6.24, delta=d)
        self.assertAlmostEqual(ps[7].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[7].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[8].a, 7.24, delta=d)
        self.assertAlmostEqual(ps[8].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[8].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[9].a, 8.24, delta=d)
        self.assertAlmostEqual(ps[9].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[9].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[10].a, 9.24, delta=d)
        self.assertAlmostEqual(ps[10].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[10].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[11].a, 10.24, delta=d)
        self.assertAlmostEqual(ps[11].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[11].T, 6.32, delta=d)
    
    def test_inclined_circular(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.00,inc=0.14,Omega=0.64) # dummy, make sure jacobi com is calculated right for higher indices
        sim.add(m=1.e-4,a=1.24,inc=0.14,Omega=0.64,l=0.632)
        sim.add(m=1.e-4,a=2.24,inc=0.14,Omega=0.64,theta=0.632)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0., delta=d)
            self.assertAlmostEqual(p.inc, 0.14, delta=d)
            self.assertAlmostEqual(p.Omega, 0.64, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[2].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[2].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[3].theta, 0.632, delta=d)
    
    def test_planar_circular(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.00) # dummy, make sure jacobi com is calculated right for higher indices
        sim.add(m=1.e-4,a=1.24,l=0.632)
        sim.add(m=1.e-4,a=2.24,theta=0.632)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0., delta=d)
            self.assertAlmostEqual(p.inc, 0., delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[2].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[2].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[3].theta, 0.632, delta=d)

    def test_retrograde_inclined_eccentric(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.24,e=0.123,inc=2.14,omega=0.12,Omega=0.64,l=0.632)
        sim.add(m=1.e-4,a=2.24,e=0.123,inc=2.14,omega=0.12,Omega=0.64,f=0.632)
        sim.add(m=1.e-4,a=3.24,e=0.123,inc=2.14,omega=0.12,Omega=0.64,theta=0.632)
        sim.add(m=1.e-4,a=4.24,e=0.123,inc=2.14,omega=0.12,Omega=0.64,M=0.632)
        sim.add(m=1.e-4,a=5.24,e=0.123,inc=2.14,omega=0.12,Omega=0.64,T=0.632)
        sim.add(m=1.e-4,a=6.24,e=0.123,inc=2.14,pomega=0.12,Omega=0.64,l=0.632)
        sim.add(m=1.e-4,a=7.24,e=0.123,inc=2.14,pomega=0.12,Omega=0.64,f=0.632)
        sim.add(m=1.e-4,a=8.24,e=0.123,inc=2.14,pomega=0.12,Omega=0.64,theta=0.632)
        sim.add(m=1.e-4,a=9.24,e=0.123,inc=2.14,pomega=0.12,Omega=0.64,M=0.632)
        sim.add(m=1.e-4,a=10.24,e=0.123,inc=2.14,pomega=0.12,Omega=0.64,T=0.632)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0.123, delta=d)
            self.assertAlmostEqual(p.inc, 2.14, delta=d)
            self.assertAlmostEqual(p.Omega, 0.64, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[1].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[1].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[1].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[2].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[2].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[2].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 3.24, delta=d)
        self.assertAlmostEqual(ps[3].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[3].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[4].a, 4.24, delta=d)
        self.assertAlmostEqual(ps[4].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[4].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[5].a, 5.24, delta=d)
        self.assertAlmostEqual(ps[5].omega, 0.12, delta=d)
        self.assertAlmostEqual(ps[5].T, 0.632, delta=d)
        self.assertAlmostEqual(ps[6].a, 6.24, delta=d)
        self.assertAlmostEqual(ps[6].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[6].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[7].a, 7.24, delta=d)
        self.assertAlmostEqual(ps[7].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[7].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[8].a, 8.24, delta=d)
        self.assertAlmostEqual(ps[8].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[8].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[9].a, 9.24, delta=d)
        self.assertAlmostEqual(ps[9].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[9].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[10].a, 10.24, delta=d)
        self.assertAlmostEqual(ps[10].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[10].T, 0.632, delta=d)
    
    def test_retrorade_planar_eccentric(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.00,e=0.123,inc=math.pi) # dummy, make sure jacobi com is calculated right for higher indices
        sim.add(m=1.e-4,a=1.24,e=0.123,inc=math.pi,pomega=0.12,l=0.632)
        sim.add(m=1.e-4,a=2.24,e=0.123,inc=math.pi,pomega=0.12,f=0.632)
        sim.add(m=1.e-4,a=3.24,e=0.123,inc=math.pi,pomega=0.12,theta=0.632)
        sim.add(m=1.e-4,a=4.24,e=0.123,inc=math.pi,pomega=0.12,M=0.632)
        sim.add(m=1.e-4,a=5.24,e=0.123,inc=math.pi,pomega=0.12,T=0.632)
        sim.add(m=1.e-4,a=6.24,e=0.123,inc=math.pi,Omega=1.25,pomega=0.12,l=0.632)
        sim.add(m=1.e-4,a=7.24,e=0.123,inc=math.pi,Omega=1.25,pomega=0.12,f=0.632)
        sim.add(m=1.e-4,a=8.24,e=0.123,inc=math.pi,Omega=1.25,pomega=0.12,theta=0.632)
        sim.add(m=1.e-4,a=9.24,e=0.123,inc=math.pi,Omega=1.25,pomega=0.12,M=0.632)
        sim.add(m=1.e-4,a=10.24,e=0.123,inc=math.pi,Omega=1.25,pomega=0.12,T=6.32)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0.123, delta=d)
            self.assertAlmostEqual(p.inc, math.pi, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[2].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[2].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[2].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[3].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[3].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[4].a, 3.24, delta=d)
        self.assertAlmostEqual(ps[4].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[4].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[5].a, 4.24, delta=d)
        self.assertAlmostEqual(ps[5].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[5].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[6].a, 5.24, delta=d)
        self.assertAlmostEqual(ps[6].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[6].T, 0.632, delta=d)
        self.assertAlmostEqual(ps[7].a, 6.24, delta=d)
        self.assertAlmostEqual(ps[7].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[7].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[8].a, 7.24, delta=d)
        self.assertAlmostEqual(ps[8].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[8].f, 0.632, delta=d)
        self.assertAlmostEqual(ps[9].a, 8.24, delta=d)
        self.assertAlmostEqual(ps[9].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[9].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[10].a, 9.24, delta=d)
        self.assertAlmostEqual(ps[10].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[10].M, 0.632, delta=d)
        self.assertAlmostEqual(ps[11].a, 10.24, delta=d)
        self.assertAlmostEqual(ps[11].pomega, 0.12, delta=d)
        self.assertAlmostEqual(ps[11].T, 6.32, delta=d)
    
    def test_retrograde_inclined_circular(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,a=1.00,inc=2.14,Omega=0.64) # dummy, make sure jacobi com is calculated right for higher indices
        sim.add(m=1.e-4,a=1.24,inc=2.14,Omega=0.64,l=0.632)
        sim.add(m=1.e-4,a=2.24,inc=2.14,Omega=0.64,theta=0.632)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0., delta=d)
            self.assertAlmostEqual(p.inc, 2.14, delta=d)
            self.assertAlmostEqual(p.Omega, 0.64, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[2].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[2].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[3].theta, 0.632, delta=d)
    
    def test_retrograde_planar_circular(self):
        sim = rebound.Simulation()
        d = 1.e-12 # abs error tolerance
        sim.add(m=1.2354)
        sim.add(m=1.e-4,inc=math.pi,a=1.00) # dummy, make sure jacobi com is calculated right for higher indices
        sim.add(m=1.e-4,inc=math.pi,a=1.24,l=0.632)
        sim.add(m=1.e-4,inc=math.pi,a=2.24,theta=0.632)

        def test_p(p):
            self.assertAlmostEqual(p.e, 0., delta=d)
            self.assertAlmostEqual(p.inc, math.pi, delta=d)

        ps = sim.particles
        for p in ps[1:]: 
            test_p(p)
           
        self.assertAlmostEqual(ps[2].a, 1.24, delta=d)
        self.assertAlmostEqual(ps[2].l, 0.632, delta=d)
        self.assertAlmostEqual(ps[2].inc, math.pi, delta=d)
        self.assertAlmostEqual(ps[3].a, 2.24, delta=d)
        self.assertAlmostEqual(ps[3].theta, 0.632, delta=d)
        self.assertAlmostEqual(ps[3].inc, math.pi, delta=d)

    def test_pal_back_and_forth(self):
        def test_p(x, y, z, vx, vy, vz):
            sim = rebound.Simulation()
            sim.add(m=1.234)
            sim.add(m=0.0002345, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
            o = sim.orbits()[0]
            sim.add(primary=sim.particles[0], m=sim.particles[1].m, a=o.a, l=o.l, h=o.pal_h, k=o.pal_k, ix=o.pal_ix, iy=o.pal_iy)
            
            d = 2.e-14 # abs error tolerance
            self.assertAlmostEqual(sim.particles[1].x, sim.particles[2].x, delta=d)
            self.assertAlmostEqual(sim.particles[1].y, sim.particles[2].y, delta=d)
            self.assertAlmostEqual(sim.particles[1].z, sim.particles[2].z, delta=d)
            self.assertAlmostEqual(sim.particles[1].vx, sim.particles[2].vx, delta=d)
            self.assertAlmostEqual(sim.particles[1].vy, sim.particles[2].vy, delta=d)
            self.assertAlmostEqual(sim.particles[1].vz, sim.particles[2].vz, delta=d)

        test_p(x=0.98, y=0.023, z=0.01, vx=-0.0151, vy=0.9981, vz=-0.01)
        test_p(x=0.198, y=0.023, z=1.01, vx=-0.01, vy=0.09981, vz=-0.01)

if __name__ == "__main__":
    unittest.main()
