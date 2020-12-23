import rebound
import unittest
import os
import rebound.data as data
from datetime import datetime

class TestMercurius(unittest.TestCase):
    
    def test_outer_solar(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        
        sim.integrator = "mercurius"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.calculate_energy()
        sim.integrate(1000)
        dE = abs((sim.calculate_energy() - E0)/E0)
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

        sim.integrate(1000)
        
        o1 = sim.particles[1].calculate_orbit(primary=sim.particles[0])
        o2 = sim.particles[2].calculate_orbit(primary=sim.particles[0])
        

        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=5e-3,a=1.0,e=0.1,primary=sim.particles[0])
        sim.add(m=1e-3,a=1.1,e=0.1,primary=sim.particles[0])

        sim.integrator = "mercurius"
        sim.dt = 1e-2*2.*3.14
        sim.N_active = 1
        sim.testparticle_type = 0

        sim.integrate(1000)
        
        o3 = sim.particles[1].calculate_orbit(primary=sim.particles[0])
        o4 = sim.particles[2].calculate_orbit(primary=sim.particles[0])

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
        
        o1 = sim.particles[1].calculate_orbit(primary=sim.particles[0])
        o2 = sim.particles[2].calculate_orbit(primary=sim.particles[0])
        

        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=5e-3,a=1.0,e=0.1,primary=sim.particles[0])
        sim.add(m=1e-3,a=1.1,e=0.1,primary=sim.particles[0])

        sim.integrator = "mercurius"
        sim.dt = 1e-2*2.*3.14
        sim.N_active = 1
        sim.testparticle_type = 1

        sim.integrate(1000)
        
        o3 = sim.particles[1].calculate_orbit(primary=sim.particles[0])
        o4 = sim.particles[2].calculate_orbit(primary=sim.particles[0])

        self.assertLess(abs(o1.e-o4.e),7e-14)
        self.assertLess(abs(o2.e-o3.e),9e-14)
        self.assertLess(abs(o1.a-o4.a),8e-14)
        self.assertLess(abs(o2.a-o3.a),4e-14)
    
    def test_outer_solar_massive(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for i in range(1,sim.N):
            sim.particles[i].m *=50.
        
        sim.integrator = "mercurius"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.calculate_energy()
        sim.integrate(1000)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,7e-8)
    
    def test_simple_collision(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        mtot0= sum([p.m for p in sim.particles])
        com0 = sim.calculate_com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.track_energy_offset = 1;
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        com1 = sim.calculate_com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=1e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=1e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=1e-16)
        mtot1= sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        dE = abs((sim.calculate_energy() - E0)/E0)
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
        sim.track_energy_offset = 1;
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        mtot1= sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        dE = abs((sim.calculate_energy() - E0)/E0)
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
        sim.track_energy_offset = 1;
        
        sim.boundary = "open"
        boxsize = 3.
        sim.configure_box(boxsize)
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,4e-6)
    
    def test_collision_with_star_simple(self):
        sim = rebound.Simulation()
        sim.add(m=1.,r=1.)
        sim.add(m=1e-3,r=1.e-3,a=0.5) 
        mtot0 = sum([p.m for p in sim.particles])
        com0 = sim.calculate_com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.track_energy_offset = 1;
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        com1 = sim.calculate_com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=1e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=1e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=1e-16)
        mtot1 = sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        self.assertEqual(N0-1,sim.N)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,1e-16)

    def test_collision_with_star(self):
        sim = rebound.Simulation()
        sim.add(m=1.,r=0.00465)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1,f=2.3) 
        sim.add(m=1e-4,r=1.4e-3,x=1.,vx=-0.4) # falling onto the star
        sim.add(m=1e-5,r=1.6e-4,a=1.5,e=0.1) 
        mtot0 = sum([p.m for p in sim.particles])
        com0 = sim.calculate_com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.track_energy_offset = 1;
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        com1 = sim.calculate_com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=1e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=1e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=1e-16)
        mtot1 = sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        self.assertEqual(N0-1,sim.N)
        dE = abs((sim.calculate_energy() - E0)/E0)
        # bad energy conservation due to democratic heliocentric!
        self.assertLess(dE,3e-2)
    
    def test_many_encounters(self):
        def get_sim():
            sim = rebound.Simulation()
            sim.add(m=1)
            for i in range(6):
                sim.add(a=1+0.2*i,e=0.1+0.1*i,f=80.*i,omega=30.*i*i,m=0.0001)
            sim.move_to_com()
            sim.dt = 0.034
            return sim
        
        sim = get_sim()
        sim.integrator = "mercurius"
        E0 = sim.calculate_energy()
        start=datetime.now()
        sim.integrate(2000)
        time_mercurius = (datetime.now()-start).total_seconds()
        dE_mercurius = abs((sim.calculate_energy() - E0)/E0)
        
        sim = get_sim()
        sim.integrator = "ias15"
        start=datetime.now()
        sim.integrate(2000)
        time_ias15 = (datetime.now()-start).total_seconds()
        dE_ias15 = abs((sim.calculate_energy() - E0)/E0)
        
        sim = get_sim()
        sim.integrator = "whfast"
        start=datetime.now()
        sim.integrate(2000)
        time_whfast = (datetime.now()-start).total_seconds()
        dE_whfast = abs((sim.calculate_energy() - E0)/E0)
        
        # Note: precision might vary on machine as initializations use cos/sin 
        # and are therefore machine dependent. 
        self.assertLess(dE_mercurius,4e-6)              # reasonable precision for mercurius
        self.assertLess(dE_mercurius/dE_whfast,1e-4)    # at least 1e4 times better than whfast
        is_travis = 'TRAVIS' in os.environ
        if not is_travis: # timing not reliable on TRAVIS
            self.assertLess(2.*time_mercurius,time_ias15) # at least 2 times faster than ias15
        


if __name__ == "__main__":
    unittest.main()

