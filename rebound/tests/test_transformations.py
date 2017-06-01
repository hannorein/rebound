import rebound
import unittest
import ctypes


def getc(sim):
    c = []
    for i in range(sim.N):
        c.append(sim.particles[0].x)
        c.append(sim.particles[0].y)
        c.append(sim.particles[0].z)
        c.append(sim.particles[0].vx)
        c.append(sim.particles[0].vy)
        c.append(sim.particles[0].vz)
        c.append(sim.particles[0].m)
        c.append(sim.particles[0].r)
    return c

class TestTransformations(unittest.TestCase):
    
    def test_democratichelio(self):
        sim = rebound.Simulation()
        sim.add(m=1.2354)
        sim.add(m=0.1,a=1.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,l=0.632)
        sim.add(m=0.01,a=5.24,e=0.2123,inc=0.014,omega=0.012,Omega=0.0164,l=10.18632)
        sim.add(m=1e-7,a=7.24,e=0.22123,inc=0.3014,omega=0.4012,Omega=0.110164,l=2.18632)
        
        elems = (rebound.Particle * sim.N)()
        p = ctypes.cast(elems,ctypes.POINTER(rebound.Particle))

        c0 = getc(sim)

        cl = rebound.clibrebound
        cl.reb_transformations_inertial_to_democraticheliocentric_posvel(sim._particles,p,sim.N)

        for i in range(sim.N):
            sim.particles[i].x = 1234.
            sim.particles[i].vx = 1234.
        cl.reb_transformations_democraticheliocentric_to_inertial_posvel(sim._particles,p,sim.N)
        
        c1 = getc(sim)
        
        for i in range(len(c0)):
            self.assertAlmostEqual(c0[i],c1[i],delta=1e-16)
        
        for i in range(sim.N):
            sim.particles[i].x = 1234.
        cl.reb_transformations_democraticheliocentric_to_inertial_pos(sim._particles,p,sim.N)
        
        c1 = getc(sim)
        
        for i in range(len(c0)):
            self.assertAlmostEqual(c0[i],c1[i],delta=1e-16)
    
    
    def test_whds(self):
        sim = rebound.Simulation()
        sim.add(m=1.2354)
        sim.add(m=0.1,a=1.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,l=0.632)
        sim.add(m=0.01,a=5.24,e=0.2123,inc=0.014,omega=0.012,Omega=0.0164,l=10.18632)
        sim.add(m=1e-7,a=7.24,e=0.22123,inc=0.3014,omega=0.4012,Omega=0.110164,l=2.18632)
        
        elems = (rebound.Particle * sim.N)()
        p = ctypes.cast(elems,ctypes.POINTER(rebound.Particle))

        c0 = getc(sim)

        cl = rebound.clibrebound
        cl.reb_transformations_inertial_to_whds_posvel(sim._particles,p,sim.N)

        for i in range(sim.N):
            sim.particles[i].x = 1234.
            sim.particles[i].vx = 1234.
        cl.reb_transformations_whds_to_inertial_posvel(sim._particles,p,sim.N)
        
        c1 = getc(sim)
        
        for i in range(len(c0)):
            self.assertAlmostEqual(c0[i],c1[i],delta=1e-16)
        
        for i in range(sim.N):
            sim.particles[i].x = 1234.
        cl.reb_transformations_whds_to_inertial_pos(sim._particles,p,sim.N)
        
        c1 = getc(sim)
        
        for i in range(len(c0)):
            self.assertAlmostEqual(c0[i],c1[i],delta=1e-16)
    

    def test_jacoobi(self):
        sim = rebound.Simulation()
        sim.add(m=1.2354)
        sim.add(m=0.1,a=1.24,e=0.123,inc=0.14,omega=0.12,Omega=0.64,l=0.632)
        sim.add(m=0.01,a=5.24,e=0.2123,inc=0.014,omega=0.012,Omega=0.0164,l=10.18632)
        sim.add(m=1e-7,a=7.24,e=0.22123,inc=0.3014,omega=0.4012,Omega=0.110164,l=2.18632)
        
        elems = (rebound.Particle * sim.N)()
        p = ctypes.cast(elems,ctypes.POINTER(rebound.Particle))
        
        elemse = (ctypes.c_double * sim.N)()

        c0 = getc(sim)

        cl = rebound.clibrebound
        
        
        cl.reb_transformations_inertial_to_jacobi_posvel(sim._particles,p,sim._particles,sim.N)

        for i in range(sim.N):
            sim.particles[i].x = 1234.
            sim.particles[i].vx = 1234.
        cl.reb_transformations_jacobi_to_inertial_posvel(sim._particles,p,sim._particles,sim.N)
        
        c1 = getc(sim)
        
        for i in range(len(c0)):
            self.assertAlmostEqual(c0[i],c1[i],delta=1e-16)
        
        for i in range(sim.N):
            sim.particles[i].x = 1234.
        cl.reb_transformations_jacobi_to_inertial_pos(sim._particles,p,sim._particles,sim.N)
        
        c1 = getc(sim)
        
        for i in range(len(c0)):
            self.assertAlmostEqual(c0[i],c1[i],delta=1e-16)
    
    
if __name__ == "__main__":
    unittest.main()
