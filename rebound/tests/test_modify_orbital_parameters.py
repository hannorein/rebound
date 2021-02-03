import rebound
import numpy as np
import unittest
import ctypes

d = 1e-14 # precision. some orbit transformations are not well behaved

class TestOrbitalElements(unittest.TestCase):
    def test_P(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].P += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].P, os0[0].P+0.5, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)
    
    def test_a(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].a += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a+0.5, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)
    
    def test_e(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].e += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e+0.5, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)
    
    def test_inc(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].inc += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc+0.5, delta=d)
    
    def test_omega(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].omega += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega+0.5, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)

    def test_pomega(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].pomega += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].pomega, os0[0].pomega+0.5, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)
   
    def test_Omega(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].Omega += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega+0.5, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)
    
    def test_f(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].f += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].f, os0[0].f+0.5, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)
    
    def test_M(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].M += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].M, os0[0].M+0.5, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)


    def test_l(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].l += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].l, os0[0].l+0.5, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)

    def test_theta(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].theta += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].theta, os0[0].theta+0.5, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)

    def test_T(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.2,Omega=0.3,inc=0.4,f=0.5)
        sim.add(m=1e-3,a=2)
        os0 = sim.calculate_orbits()
        sim.particles[1].T += 0.5
        os1 = sim.calculate_orbits()
        self.assertAlmostEqual(os1[0].a, os0[0].a, delta=d)
        self.assertAlmostEqual(os1[0].e, os0[0].e, delta=d)
        self.assertAlmostEqual(os1[0].T, os0[0].T+0.5, delta=d)
        self.assertAlmostEqual(os1[0].omega, os0[0].omega, delta=d)
        self.assertAlmostEqual(os1[0].Omega, os0[0].Omega, delta=d)
        self.assertAlmostEqual(os1[0].inc, os0[0].inc, delta=d)



if __name__ == "__main__":
    unittest.main()
