import rebound
import unittest
import sys
import warnings
import os
from datetime import datetime

class TestMercurius(unittest.TestCase):
    
    def test_no_effect_tp(self):
        # tests if test particle encounters have an effect
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.move_to_com()
        sim.N_active=2
        sim.integrator = "mercurius"
        sim.dt = 0.1
        sim2 = sim.copy()
        p = sim.particles[1].copy()
        p.x += 0.01
        p.m = 0
        sim.add(p)
        sim.step()
        sim2.step()

        self.assertEqual(sim.particles[1].x,sim2.particles[1].x)
        self.assertEqual(sim.particles[1].vx,sim2.particles[1].vx)
        self.assertEqual(sim.particles[0].x,sim2.particles[0].x)
        self.assertEqual(sim.particles[0].vx,sim2.particles[0].vx)

    def test_outer_solar(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        
        sim.integrator = "mercurius"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.energy()
        sim.integrate(1000)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,2e-10)
    
    def test_order_doesnt_matter_tp0(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1.1,e=0.1,primary=sim.particles[0])
        sim.add(m=5e-3,a=1.0,e=0.1,primary=sim.particles[0])

        sim.integrator = "mercurius"
        #sim.integrator = "whfast"
        #sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 1e-2*2.*3.14
        sim.N_active = 1
        sim.testparticle_type = 0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.integrate(1000)
            self.assertEqual(1,len(w))
        
        o1 = sim.particles[1].orbit(primary=sim.particles[0])
        o2 = sim.particles[2].orbit(primary=sim.particles[0])
        

        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=5e-3,a=1.0,e=0.1,primary=sim.particles[0])
        sim.add(m=1e-3,a=1.1,e=0.1,primary=sim.particles[0])

        sim.integrator = "mercurius"
        sim.dt = 1e-2*2.*3.14
        sim.N_active = 1
        sim.testparticle_type = 0

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.integrate(1000)
            self.assertEqual(1,len(w))
        
        o3 = sim.particles[1].orbit(primary=sim.particles[0])
        o4 = sim.particles[2].orbit(primary=sim.particles[0])

        self.assertLess(abs(o1.e-o4.e),2e-16)
        self.assertLess(abs(o2.e-o3.e),2e-16)
        self.assertLess(abs(o1.a-o4.a),2e-16)
        self.assertLess(abs(o2.a-o3.a),2e-16)
    
    def test_order_doesnt_matter_tp1(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1.1,e=0.1,primary=sim.particles[0])
        sim.add(m=5e-3,a=1.0,e=0.1,primary=sim.particles[0])

        sim.integrator = "mercurius"
        sim.dt = 1e-2*2.*3.14
        sim.N_active = 1
        sim.testparticle_type = 1

        sim.integrate(1000)
        
        o1 = sim.particles[1].orbit(primary=sim.particles[0])
        o2 = sim.particles[2].orbit(primary=sim.particles[0])
        

        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=5e-3,a=1.0,e=0.1,primary=sim.particles[0])
        sim.add(m=1e-3,a=1.1,e=0.1,primary=sim.particles[0])

        sim.integrator = "mercurius"
        sim.dt = 1e-2*2.*3.14
        sim.N_active = 1
        sim.testparticle_type = 1

        sim.integrate(1000)
        
        o3 = sim.particles[1].orbit(primary=sim.particles[0])
        o4 = sim.particles[2].orbit(primary=sim.particles[0])

        self.assertLess(abs(o1.e-o4.e),2e-13)
        self.assertLess(abs(o2.e-o3.e),9e-14)
        self.assertLess(abs(o1.a-o4.a),9e-14)
        self.assertLess(abs(o2.a-o3.a),4e-14)
    
    def test_outer_solar_massive(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for i in range(1,sim.N):
            sim.particles[i].m *=50.
        
        sim.integrator = "mercurius"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.energy()
        sim.integrate(1000)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,7e-8)
    
    def test_simple_collision(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        mtot0= sum([p.m for p in sim.particles])
        com0 = sim.com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.track_energy_offset = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.energy()
        sim.integrate(1)
        com1 = sim.com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=1e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=1e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=1e-16)
        mtot1= sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,3e-9)
        self.assertEqual(N0-1,sim.N)


    def test_planetesimal_collision(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.N_active = 2
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        mtot0= sum([p.m for p in sim.particles])
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.testparticle_type = 1
        sim.track_energy_offset = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.energy()
        sim.integrate(1)
        mtot1= sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,3e-9)
        self.assertEqual(N0-1,sim.N)

    def test_massive_ejection(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-4,r=1.6e-4,a=0.5,e=0.1)
        sim.add(m=1e-6,r=4e-5,a=0.6)
        sim.particles[2].vy *= 2

        sim.N_active = 2
        
        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.testparticle_type = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.track_energy_offset = 1
        
        sim.boundary = "open"
        boxsize = 3.
        sim.configure_box(boxsize)
        
        E0 = sim.energy()
        sim.integrate(1)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,4e-6)
    
    def test_collision_with_star_simple(self):
        sim = rebound.Simulation()
        sim.add(m=1.,r=1.)
        sim.add(m=1e-3,r=1.e-3,a=0.5) 
        mtot0 = sum([p.m for p in sim.particles])
        com0 = sim.com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.track_energy_offset = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.energy()
        sim.integrate(1)
        com1 = sim.com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=1e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=1e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=1e-16)
        mtot1 = sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        self.assertEqual(N0-1,sim.N)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,1e-16)

    def test_collision_with_star(self):
        sim = rebound.Simulation()
        sim.add(m=1.,r=0.00465)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1,f=2.3) 
        sim.add(m=1e-4,r=1.4e-3,x=1.,vx=-0.4) # falling onto the star
        sim.add(m=1e-5,r=1.6e-4,a=1.5,e=0.1) 
        mtot0 = sum([p.m for p in sim.particles])
        com0 = sim.com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.track_energy_offset = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.energy()
        sim.integrate(1)
        com1 = sim.com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=1e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=1e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=1e-16)
        mtot1 = sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        self.assertEqual(N0-1,sim.N)
        dE = abs((sim.energy() - E0)/E0)
        # bad energy conservation due to democratic heliocentric!
        self.assertLess(dE,3e-2)
    
    def test_many_encounters(self):
        def get_sim():
            sim = rebound.Simulation()
            sim.add(m=1)
            # Setup using xyz instead of orbital elements for
            # machine independent test
            sim.add(m=0.0001,x=0.90000, y=0.00000, vx=0.00000, vy=1.10360)
            sim.add(m=0.0001, x=-1.17676, y=-0.05212, vx=0.22535, vy=-0.90102)
            sim.add(m=0.0001, x=-1.66025, y=-0.69852, vx=0.18932, vy=-0.60030)
            sim.add(m=0.0001, x=0.57904, y=1.03836, vx=-0.69267, vy=0.75995)
            sim.add(m=0.0001, x=-0.41683, y=0.83128, vx=-1.03478, vy=-0.72482)
            sim.add(m=0.0001, x=1.83969, y=0.32938, vx=-0.55114, vy=0.51646)
            sim.move_to_com()
            sim.dt = 0.034
            return sim
        
        sim = get_sim()
        sim.integrator = "mercurius"
        E0 = sim.energy()
        start=datetime.now()
        sim.integrate(2000)
        time_mercurius = (datetime.now()-start).total_seconds()
        dE_mercurius = abs((sim.energy() - E0)/E0)
        
        sim = get_sim()
        sim.integrator = "ias15"
        start=datetime.now()
        sim.integrate(2000)
        time_ias15 = (datetime.now()-start).total_seconds()
        dE_ias15 = abs((sim.energy() - E0)/E0)
        
        sim = get_sim()
        sim.integrator = "whfast"
        start=datetime.now()
        sim.integrate(2000)
        time_whfast = (datetime.now()-start).total_seconds()
        dE_whfast = abs((sim.energy() - E0)/E0)
        
        # Note: precision might vary on machine as initializations use cos/sin 
        # and are therefore machine dependent. 
        self.assertLess(dE_mercurius,5e-6)              # reasonable precision for mercurius
        self.assertLess(dE_mercurius/dE_whfast,1e-4)    # at least 1e4 times better than whfast
        if os.getenv("CI") != "true":
            self.assertLess(time_mercurius,5.0*time_ias15)  # not much slower than ias15 (often fails in unit tests)
        if sys.maxsize > 2**32: # 64 bit
            self.assertEqual(7060.644251181158, sim.particles[5].x) # Check if bitwise unchanged
        


if __name__ == "__main__":
    unittest.main()

