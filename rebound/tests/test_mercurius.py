import rebound
import unittest
import os
import rebound.data as data

class TestMercurius(unittest.TestCase):
    
    def test_outer_solar(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        
        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.calculate_energy()
        sim.integrate(1000)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,2e-10)
    
    def test_outer_solar_massive(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for i in range(1,sim.N):
            sim.particles[i].m *=50.
        
        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.calculate_energy()
        sim.integrate(1000)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,7e-8)
        self.assertGreater(sim.ri_mercurius._encounterAllocatedN,1)

    def test_planetesimal_collision(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.N_active = 2
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        sim.move_to_com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 0.0001
        sim.testparticle_type = 1
        sim.track_energy_offset = 1;
        sim.collision_resolve_keep_sorted = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,2e-9)
        self.assertEqual(N0-1,sim.N)

    def test_massive_ejection(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-4,r=1.6e-4,a=0.5,e=0.1)
        sim.add(m=1e-6,r=4e-5,a=0.6)
        sim.particles[2].vy *= 2
        sim.N_active = 2
        sim.move_to_com()
        
        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 0.001
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
        self.assertLess(dE,5e-11)


class TestMercuriusWHDS(unittest.TestCase):
    
    def test_outer_solar(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        
        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "whds"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.calculate_energy()
        sim.integrate(1000)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,5e-10)
    
    def test_outer_solar_massive(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for i in range(1,sim.N):
            sim.particles[i].m *=50.
        
        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "whds"
        P = sim.particles[1].P
        sim.dt = 1e-3*P
        
        E0 = sim.calculate_energy()
        sim.integrate(1000)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,8e-8)
        self.assertGreater(sim.ri_mercurius._encounterAllocatedN,1)

    def test_planetesimal_collision(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.N_active = 2
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        sim.move_to_com()
        N0 = sim.N

        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "whds"
        sim.dt = 0.0001
        sim.testparticle_type = 1
        sim.track_energy_offset = 1;
        sim.collision_resolve_keep_sorted = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,5e-8)
        self.assertEqual(N0-1,sim.N)

    def test_massive_ejection(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-4,r=1.6e-4,a=0.5,e=0.1)
        sim.add(m=1e-6,r=4e-5,a=0.6)
        sim.particles[2].vy *= 2
        sim.N_active = 2
        sim.move_to_com()
        
        sim.integrator = "mercurius"
        sim.ri_whfast.coordinates = "whds"
        sim.dt = 0.001
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
        self.assertLess(dE,5e-11)

if __name__ == "__main__":
    unittest.main()

