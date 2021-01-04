import rebound
import unittest
import math
import rebound.data
import warnings

coordinatelist = ["democraticheliocentric","whds","jacobi"]
class TestIntegratorWHFastTestParticle(unittest.TestCase):
    pass

def create_whfast_testparticle(coordinates, N, N_active):
    def do_test(self):
        sim = rebound.Simulation()
        sim.ri_whfast.coordinates = coordinates
        sim.integrator = "whfast"
        sim.dt=1e-1
        sim.add(m=1)
        for i in range(N):
            sim.add(m=0,P=1,e=0.1,f=i)
  
        sim2 = sim.copy()
        if N_active>0:
            sim2.N_active = N_active
        sim2.integrate(1)

        # some combinations can't do a simple keplerian orbit exactly (sad but true)
        eps = 1e-14
        for i in range(sim.N):
            self.assertLess(abs(sim.particles[i].x-sim2.particles[i].x),eps)
            self.assertLess(abs(sim.particles[i].vx-sim2.particles[i].vx),eps)
            self.assertLess(abs(sim.particles[i].y-sim2.particles[i].y),eps)
            self.assertLess(abs(sim.particles[i].vy-sim2.particles[i].vy),eps)
        
        # all should be able to do this exactly
        sim.integrate(1)
        eps = 1e-16
        for i in range(sim.N):
            self.assertLess(abs(sim.particles[i].x-sim2.particles[i].x),eps)
            self.assertLess(abs(sim.particles[i].vx-sim2.particles[i].vx),eps)
            self.assertLess(abs(sim.particles[i].y-sim2.particles[i].y),eps)
            self.assertLess(abs(sim.particles[i].vy-sim2.particles[i].vy),eps)
    return do_test


def create_whfast_testparticle_withplanet(coordinates, N, N_active):
    def do_test(self):
        eps = 1e-13
        sim = rebound.Simulation()
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = coordinates
        sim.dt=1e-3
        sim.add(m=1)
        sim.add(m=1e-3,P=0.4)
        sim.add(m=1e-3,P=0.7)
        for i in range(N):
            sim.add(m=0,P=1,e=0.1,f=i)
  
        sim2 = sim.copy()
        if N_active>0:
            sim2.N_active = N_active
       
        sim2.integrate(1)
        sim.integrate(1)
        for i in range(sim.N):
            self.assertLess(abs(sim.particles[i].x-sim2.particles[i].x),eps)
            self.assertLess(abs(sim.particles[i].vx-sim2.particles[i].vx),eps)
            self.assertLess(abs(sim.particles[i].y-sim2.particles[i].y),eps)
            self.assertLess(abs(sim.particles[i].vy-sim2.particles[i].vy),eps)
    return do_test

def create_whfast_testparticletype1(coordinates, N_active):
    def do_test(self):
        eps = 1e-16
        sim = rebound.Simulation()
        sim.ri_whfast.coordinates = coordinates
        sim.testparticle_type = 1
        sim.integrator = "whfast"
        sim.dt=1e-3
        sim.add(m=1)
        sim.add(m=1e-3,P=1,e=0.1)
  
        sim2 = sim.copy()
        if N_active>0:
            sim2.N_active = N_active
        sim2.integrate(1)
        sim.integrate(1)

        for i in range(sim.N):
            self.assertLess(abs(sim.particles[i].x-sim2.particles[i].x),eps)
            self.assertLess(abs(sim.particles[i].vx-sim2.particles[i].vx),eps)
            self.assertLess(abs(sim.particles[i].y-sim2.particles[i].y),eps)
            self.assertLess(abs(sim.particles[i].vy-sim2.particles[i].vy),eps)
    return do_test
def create_whfast_testparticletype1_withplanet(coordinates, N_active):
    def do_test(self):
        eps = 1e-16
        sim = rebound.Simulation()
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = coordinates
        sim.testparticle_type = 1
        sim.dt=1e-3
        sim.add(m=1)
        sim.add(m=1e-3,P=0.4)
        sim.add(m=1e-3,P=0.7,e=0.1)
        sim.add(m=1e-3,P=1.0)
  
        sim2 = sim.copy()
        if N_active>0:
            sim2.N_active = N_active
       
        sim2.integrate(1)
        sim.integrate(1)
        for i in range(sim.N):
            self.assertLess(abs(sim.particles[i].x-sim2.particles[i].x),eps)
            self.assertLess(abs(sim.particles[i].vx-sim2.particles[i].vx),eps)
            self.assertLess(abs(sim.particles[i].y-sim2.particles[i].y),eps)
            self.assertLess(abs(sim.particles[i].vy-sim2.particles[i].vy),eps)
    return do_test

## Testparticles with mass currently lead to unexpexted behaviour:
def create_whfast_massivetestparticle(coordinates, N):
    def do_test(self):
        eps = 2e-13
        sim = rebound.Simulation()
        sim.ri_whfast.coordinates = coordinates
        sim.integrator = "whfast"
        sim.dt=1e-3
        sim.add(m=1)
        sim.add(m=1e-3,P=0.4)
        for i in range(N):
            sim.add(m=0,P=1,e=0.1,f=i)
  
        sim2 = sim.copy()
        for i in range(N):
            sim2.particles[i+2].m = 1 # particles have zero mass for sim, but finite for sim2

        sim2.N_active = 2

        sim.integrate(1)
        sim2.integrate(1)

        for i in range(sim.N):
            self.assertLess(abs(sim.particles[i].x-sim2.particles[i].x),eps)
            self.assertLess(abs(sim.particles[i].vx-sim2.particles[i].vx),eps)
            self.assertLess(abs(sim.particles[i].y-sim2.particles[i].y),eps)
            self.assertLess(abs(sim.particles[i].vy-sim2.particles[i].vy),eps)
    return do_test


for N in [1,2]: 
    for coordinates in coordinatelist:
        for N_active in [-1]+list(range(1,N+2)):
            test_method = create_whfast_testparticle(coordinates,N, N_active)
            test_method.__name__ = "test_whfast_testparticle_N%d_Nactive%d_"%(N,N_active)+coordinates
            setattr(TestIntegratorWHFastTestParticle, test_method.__name__, test_method)
            
        for N_active in [-1]+list(range(3,N+4)):
            test_method = create_whfast_testparticle_withplanet(coordinates,N, N_active)
            test_method.__name__ = "test_whfast_testparticle_withplanet_N%d_Nactive%d_"%(N,N_active)+coordinates
            setattr(TestIntegratorWHFastTestParticle, test_method.__name__, test_method)
        test_method = create_whfast_massivetestparticle(coordinates,N)
        test_method.__name__ = "test_whfast_massivetestparticle_N%d_"%(N)+coordinates
        setattr(TestIntegratorWHFastTestParticle, test_method.__name__, test_method)

for coordinates in coordinatelist:
    for N_active in [-1,1]:
        test_method = create_whfast_testparticletype1(coordinates, N_active)
        test_method.__name__ = "test_whfast_testparticletype1_Nactive%d_"%(N_active)+coordinates
        setattr(TestIntegratorWHFastTestParticle, test_method.__name__, test_method)
    for N_active in [-1,3]:
        test_method = create_whfast_testparticletype1_withplanet(coordinates, N_active)
        test_method.__name__ = "test_whfast_testparticletype1_withplanet_Nactive%d_"%(N_active)+coordinates
        setattr(TestIntegratorWHFastTestParticle, test_method.__name__, test_method)

if __name__ == "__main__":
    unittest.main()
