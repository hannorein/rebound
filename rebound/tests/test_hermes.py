import rebound
import unittest
import os
import rebound.data as data
import numpy as np
from ctypes import byref

class TestHermes(unittest.TestCase):
    def test_no_close_encounter(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-3, a=1.523,e=0.0146,f=0.24)
        sim.add(m=1.e-3, a=2.423523,e=0.01246,f=0.324)
        sim.integrator = "hermes"
        P = sim.particles[1].P
        sim.dt = 1e-4*P
        sim.integrate(P)
        x_hermes = sim.particles[1].x
        
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-3, a=1.523,e=0.0146,f=0.24)
        sim.add(m=1.e-3, a=2.423523,e=0.01246,f=0.324)
        sim.integrator = "whfasthelio"
        P = sim.particles[1].P
        sim.dt = 1e-4*P
        sim.integrate(P)
        x_whfast = sim.particles[1].x

        self.assertEqual(x_hermes,x_whfast)

    def test_close_encounter(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-3, a=1.523,e=0.0146,f=0.24)
        rh = sim.particles[1].a*pow(sim.particles[1].m/(3.*sim.particles[0].m),1./3)
        dust = rebound.Particle(simulation=sim, primary=sim.particles[1], a=0.25*rh, e=0.000123, f=2.3, m=1e-8)
        sim.add(dust)
        sim.integrator = "hermes"
        sim.gravity = "basic"
        sim.ri_hermes.hill_switch_factor = 2.
        P = sim.particles[1].P
        sim.dt = 1e-4*P
        for i in range(1000):
            sim.step()
        x_hermes = sim.particles[1].x
        
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1.e-3, a=1.523,e=0.0146,f=0.24)
        dust = rebound.Particle(simulation=sim, primary=sim.particles[1], a=0.25*rh, e=0.000123, f=2.3, m=1e-8)
        sim.add(dust)
        sim.integrator = "ias15"
        sim.gravity = "basic"
        P = sim.particles[1].P
        dt0 = 1e-4*P
        sim.dt = dt0
        t = 0.
        for i in range(1000):
            t += dt0/2.
            t += dt0/2.
            sim.integrate(t)
        x_ias15 = sim.particles[1].x
        self.assertEqual(x_hermes,x_ias15)

    def test_planetesimal_collision(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.N_active = 2
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        sim.move_to_com()
        
        sim.integrator = "hermes"
        #sim.gravity = "basic"
        sim.ri_hermes.hill_switch_factor = 3.
        sim.ri_hermes.solar_switch_factor = 20.
        sim.ri_hermes.adaptive_hill_switch_factor = 0
        sim.dt = 0.0001
        sim.testparticle_type = 1
        sim.track_energy_offset = 1;
        sim.collision_resolve_keep_sorted = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        
        E0 = sim.calculate_energy()
        sim.integrate(1)
        dE = abs((sim.calculate_energy() - E0)/E0)
        self.assertLess(dE,5e-13)

    def test_massive_ejection(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-4,r=1.6e-4,a=0.5,e=0.1)
        sim.add(m=1e-6,r=4e-5,a=0.6)
        sim.particles[2].vy *= 2
        sim.N_active = 2
        sim.move_to_com()
        
        sim.integrator = "hermes"
        #sim.gravity = "basic"
        sim.ri_hermes.hill_switch_factor = 3.
        sim.ri_hermes.solar_switch_factor = 20.
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

    def test_autoHSF_case_total_overlap(self):
        #test 1 - both massive bodies, outer body completely overlaps inner - CASE 2
        sim = rebound.Simulation()
        sim.integrator = "hermes"
        sim.ri_hermes.hill_switch_factor = 1.
        sim.dt = 0.001
        sim.add(m=1)
        sim.add(m=1e-8, a=1.)
        sim.add(m=1e-8, a=1.05, e=0.75)
        sim.move_to_com()
        sim.integrate(3)
        self.assertGreater(sim.ri_hermes._current_hill_switch_factor,1)

        #test 2 - outer body is test - CASE 2
        sim2 = rebound.Simulation()
        sim2.integrator = "hermes"
        sim2.ri_hermes.hill_switch_factor = 1.
        sim2.testparticle_type = 1
        sim2.dt = 0.001
        sim2.add(m=1)
        sim2.add(m=1e-8, a=1.)
        sim2.N_active = sim2.N
        sim2.add(m=1e-8, a=1.1, e=0.75)
        sim2.move_to_com()
        sim2.integrate(3)
        self.assertEqual(sim2.ri_hermes._current_hill_switch_factor,sim.ri_hermes._current_hill_switch_factor)

        #test 3 - inner body now has large eccentricity, completely overlaps outer - CASE 1
        sim3 = rebound.Simulation()
        sim3.integrator = "hermes"
        sim3.ri_hermes.hill_switch_factor = 1.
        sim3.dt = 0.001
        sim3.add(m=1)
        sim3.add(m=1e-8, a=1., e=0.75)
        sim3.N_active = sim3.N
        sim3.add(m=1e-8, a=1.05)
        sim3.move_to_com()
        sim3.integrate(3)
        self.assertEqual(sim3.ri_hermes._current_hill_switch_factor,sim.ri_hermes._current_hill_switch_factor)

    def test_autoHSF_case_partial_overlap(self):
        #test 1 - both massive bodies, partial overlap between inner and outer body - CASE 3
        sim = rebound.Simulation()
        sim.integrator = "hermes"
        sim.ri_hermes.hill_switch_factor = 1.
        sim.dt = 0.001
        sim.add(m=1)
        sim.add(m=1e-8, a=1., e=0.4)
        sim.add(m=1e-8, a=1.05, e=0.4)
        sim.move_to_com()
        sim.integrate(3)
        self.assertGreater(sim.ri_hermes._current_hill_switch_factor,1)

        #test 2 - partial overlap between inner (test) and outer (massive) body - CASE 4
        sim2 = rebound.Simulation()
        sim2.integrator = "hermes"
        sim2.ri_hermes.hill_switch_factor = 1.
        sim2.testparticle_type = 1
        sim2.dt = 0.001
        sim2.add(m=1)
        sim2.add(m=1e-8, a=1.05, e=0.4)
        sim2.N_active = sim2.N
        sim2.add(m=1e-8, a=1.0, e=0.4)
        sim2.move_to_com()
        sim2.integrate(3)
        self.assertEqual(sim2.ri_hermes._current_hill_switch_factor,sim.ri_hermes._current_hill_switch_factor)
'''
    def test_autoHSF_proper_calculation(self):
        sim = rebound.Simulation()
        sim.integrator="hermes"
        sim.ri_hermes.hill_switch_factor = 1.
        sim.dt = 0.001
        sim.add(m=1)
        sim.add(a=1, m=1e-11)
        sim.add(a=1, m=1e-11, e=0.02)
        dv = np.sqrt((sim.particles[1].vx - sim.particles[2].vx)**2 + (sim.particles[1].vy - sim.particles[2].vy)**2)
        rhillsum = (sim.particles[1].a + sim.particles[2].a)*np.power(sim.particles[1].m/(3.*sim.particles[0].m),1./3.)
        
        sim.integrate(sim.dt)
        HSF_theory = 1.25**(np.ceil(np.log10(sim.dt/(0.25*rhillsum/dv)/np.log10(1.25))))
        print HSF_theory, sim.ri_hermes._current_hill_switch_factor

#        print("1/4 encounter timescale:", 0.25*2*rhillsum/dv)
#        print("timestep               :", sim.dt)
#        print("should have HSF > than :", sim.dt/(0.25*rhillsum/dv))
#        print("HSF returned by hermes :", sim.ri_hermes._current_hill_switch_factor)
        #rebound.clibrebound.reb_integrator_hermes_check_for_encounter(ctypes.byref(sim))
        #print("HSF returned by hermes :", sim.ri_hermes._current_hill_switch_factor)
'''

if __name__ == "__main__":
    unittest.main()

