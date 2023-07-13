import rebound
import unittest
import datetime
import socket
import warnings

class TestCopy(unittest.TestCase):
    def test_copy_is_same(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.integrate(42.)
        sim.save("test.bin")
        

        sim_copy = sim.copy()
        
        self.assertEqual(sim.t,sim_copy.t)
        for i in range(sim.N):
            self.assertEqual(sim.particles[i].x,sim_copy.particles[i].x)
            self.assertEqual(sim.particles[i].y,sim_copy.particles[i].y)
            self.assertEqual(sim.particles[i].z,sim_copy.particles[i].z)
            self.assertEqual(sim.particles[i].vx,sim_copy.particles[i].vx)
            self.assertEqual(sim.particles[i].vy,sim_copy.particles[i].vy)
            self.assertEqual(sim.particles[i].vz,sim_copy.particles[i].vz)
        
        sim_copy.integrate(84.)
        sim.integrate(84.)

        self.assertEqual(sim.t,sim_copy.t)
        for i in range(sim.N):
            self.assertEqual(sim.particles[i].x,sim_copy.particles[i].x)
            self.assertEqual(sim.particles[i].y,sim_copy.particles[i].y)
            self.assertEqual(sim.particles[i].z,sim_copy.particles[i].z)
            self.assertEqual(sim.particles[i].vx,sim_copy.particles[i].vx)
            self.assertEqual(sim.particles[i].vy,sim_copy.particles[i].vy)
            self.assertEqual(sim.particles[i].vz,sim_copy.particles[i].vz)
        
        sim.integrate(126.)
        
        self.assertNotEqual(sim.t,sim_copy.t)
        for i in range(sim.N):
            self.assertNotEqual(sim.particles[i].x,sim_copy.particles[i].x)
            self.assertNotEqual(sim.particles[i].y,sim_copy.particles[i].y)
            self.assertNotEqual(sim.particles[i].z,sim_copy.particles[i].z)
            self.assertNotEqual(sim.particles[i].vx,sim_copy.particles[i].vx)
            self.assertNotEqual(sim.particles[i].vy,sim_copy.particles[i].vy)
            self.assertNotEqual(sim.particles[i].vz,sim_copy.particles[i].vz)

class TestMultiply(unittest.TestCase):
    def test_multiply_with_minus_one(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim1.integrate(2.)
       
        sim2 = sim1*(-1.)
        
        for i in range(sim1.N):
            self.assertEqual(-sim1.particles[i].x,sim2.particles[i].x)
            self.assertEqual(-sim1.particles[i].y,sim2.particles[i].y)
            self.assertEqual(-sim1.particles[i].z,sim2.particles[i].z)
            self.assertEqual(-sim1.particles[i].vx,sim2.particles[i].vx)
            self.assertEqual(-sim1.particles[i].vy,sim2.particles[i].vy)
            self.assertEqual(-sim1.particles[i].vz,sim2.particles[i].vz)
    
    def test_rmultiply_with_minus_one(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim1.integrate(2.)
       
        sim2 = -1.*sim1
        
        for i in range(sim1.N):
            self.assertEqual(-sim1.particles[i].x,sim2.particles[i].x)
            self.assertEqual(-sim1.particles[i].y,sim2.particles[i].y)
            self.assertEqual(-sim1.particles[i].z,sim2.particles[i].z)
            self.assertEqual(-sim1.particles[i].vx,sim2.particles[i].vx)
            self.assertEqual(-sim1.particles[i].vy,sim2.particles[i].vy)
            self.assertEqual(-sim1.particles[i].vz,sim2.particles[i].vz)
    
    def test_multiply_with_zero(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim1.integrate(2.)
        sim2 = 0*sim1
        
        for i in range(sim1.N):
            self.assertNotEqual(-sim1.particles[i].x,0.)
            self.assertNotEqual(-sim1.particles[i].y,0.)
            self.assertNotEqual(-sim1.particles[i].z,0.)
            self.assertNotEqual(-sim1.particles[i].vx,0.)
            self.assertNotEqual(-sim1.particles[i].vy,0.)
            self.assertNotEqual(-sim1.particles[i].vz,0.)
        for i in range(sim2.N):
            self.assertEqual(-sim2.particles[i].x,0.)
            self.assertEqual(-sim2.particles[i].y,0.)
            self.assertEqual(-sim2.particles[i].z,0.)
            self.assertEqual(-sim2.particles[i].vx,0.)
            self.assertEqual(-sim2.particles[i].vy,0.)
            self.assertEqual(-sim2.particles[i].vz,0.)
    
    def test_div2(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim1.integrate(2.)
        sim2 = sim1/2
        
        for i in range(sim1.N):
            self.assertEqual(2*sim2.particles[i].x,sim1.particles[i].x)
            self.assertEqual(2*sim2.particles[i].y,sim1.particles[i].y)
            self.assertEqual(2*sim2.particles[i].z,sim1.particles[i].z)
            self.assertEqual(2*sim2.particles[i].vx,sim1.particles[i].vx)
            self.assertEqual(2*sim2.particles[i].vy,sim1.particles[i].vy)
            self.assertEqual(2*sim2.particles[i].vz,sim1.particles[i].vz)

class TestAdd(unittest.TestCase):
    def test_add(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim2 = sim1.copy()
        sim1.integrate(2.)

        sim3 = sim1 + sim2
        
        for i in range(sim3.N):
            self.assertEqual(sim3.particles[i].x, sim1.particles[i].x +sim2.particles[i].x )
            self.assertEqual(sim3.particles[i].y, sim1.particles[i].y +sim2.particles[i].y )
            self.assertEqual(sim3.particles[i].z, sim1.particles[i].z +sim2.particles[i].z )
            self.assertEqual(sim3.particles[i].vx,sim1.particles[i].vx+sim2.particles[i].vx)
            self.assertEqual(sim3.particles[i].vy,sim1.particles[i].vy+sim2.particles[i].vy)
            self.assertEqual(sim3.particles[i].vz,sim1.particles[i].vz+sim2.particles[i].vz)
    
    def test_iadd(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim2 = sim1.copy()
        sim1.integrate(2.)

        sim3 = sim1 + sim2
        sim1 += sim2
        
        for i in range(sim3.N):
            self.assertEqual(sim3.particles[i].x, sim1.particles[i].x )
            self.assertEqual(sim3.particles[i].y, sim1.particles[i].y )
            self.assertEqual(sim3.particles[i].z, sim1.particles[i].z )
            self.assertEqual(sim3.particles[i].vx,sim1.particles[i].vx)
            self.assertEqual(sim3.particles[i].vy,sim1.particles[i].vy)
            self.assertEqual(sim3.particles[i].vz,sim1.particles[i].vz)
    
class TestSubtract(unittest.TestCase):
    def test_subtract(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim2 = sim1.copy()
        sim1.integrate(2.)
        
        sim3 = sim1 - sim2
        
        for i in range(sim3.N):
            self.assertEqual(sim3.particles[i].x, sim1.particles[i].x -sim2.particles[i].x )
            self.assertEqual(sim3.particles[i].y, sim1.particles[i].y -sim2.particles[i].y )
            self.assertEqual(sim3.particles[i].z, sim1.particles[i].z -sim2.particles[i].z )
            self.assertEqual(sim3.particles[i].vx,sim1.particles[i].vx-sim2.particles[i].vx)
            self.assertEqual(sim3.particles[i].vy,sim1.particles[i].vy-sim2.particles[i].vy)
            self.assertEqual(sim3.particles[i].vz,sim1.particles[i].vz-sim2.particles[i].vz)
    
    def test_isubtract(self):
        sim1 = rebound.Simulation()
        sim1.add(m=1)
        sim1.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim1.integrator = "whfast"
        sim2 = sim1.copy()
        sim1.integrate(2.)
        
        sim3 = sim1 - sim2
        sim1 -= sim2
        
        for i in range(sim3.N):
            self.assertEqual(sim3.particles[i].x, sim1.particles[i].x )
            self.assertEqual(sim3.particles[i].y, sim1.particles[i].y )
            self.assertEqual(sim3.particles[i].z, sim1.particles[i].z )
            self.assertEqual(sim3.particles[i].vx,sim1.particles[i].vx)
            self.assertEqual(sim3.particles[i].vy,sim1.particles[i].vy)
            self.assertEqual(sim3.particles[i].vz,sim1.particles[i].vz)
    

if __name__ == "__main__":
    unittest.main()
