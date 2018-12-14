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


if __name__ == "__main__":
    unittest.main()
