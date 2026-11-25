import ctypes
import rebound
import unittest
import numpy as np
import math
import sys
import warnings
import os
from datetime import datetime

def collision_add_particle(sim_pointer, collision):
    sim = sim_pointer.contents
    sim.add(m=1e-10, a=1.) # meaningless
    sim.add(m=1e-10, a=2.) # meaningless
    sim.add(m=1e-10, a=3.) # meaningless
    return 3

def peri_switching_off(sim):
    # This reproduces reb_integrator_trace_switch_peri_none, since functionality is removed in REBOUND v5
    sim.integrator.peri_crit_eta = 1e100 # pericenter criterion never fires
    sim.integrator.r_crit_WB = 0. # companion check never fires
    return sim

def chaotic_exchange_sim():
    sim = rebound.Simulation()
    # Setup using xyz instead of orbital elements for
    # machine independent test

    star_m = 1
    jup_m = 0.01 / (star_m - 0.01)
    jup_a = 5.2
    jup_e = 0.0

    star_x = -(jup_m / (star_m + jup_m)) * (jup_a * (1 + jup_e))
    star_vy = -(jup_m / (star_m + jup_m)) * np.sqrt(((star_m + jup_m) / jup_a) * ((1 - jup_e) / (1 + jup_e)))

    jup_x = (star_m / (star_m + jup_m)) * (jup_a * (1 + jup_e))
    jup_vy = (star_m / (star_m + jup_m)) * np.sqrt(((star_m + jup_m) / jup_a) * ((1 - jup_e) / (1 + jup_e)))

    t_x = 4.42
    t_vy = 0.0072 * (365.25) * (1 / (2 * np.pi))

    sim = rebound.Simulation()
    sim.add(m=star_m, x=star_x, vy=star_vy)
    sim.add(m=jup_m, x=jup_x, vy = jup_vy)
    sim.add(m=0, x=t_x + star_x, vy = t_vy + star_vy)
    sim.add(m=1.0, x=100, vy=0.1, name="widebinary")
    sim.integrator = "trace"
    sim.integrator.coordinates = "widebinary"
    sim.integrator.r_crit_hill = 5
    sim.dt = (8./365.)*2.*math.pi
    return sim

def jupsat_sim(wb_companion=True, companion_first=False):
    sim = rebound.Simulation()
    sun = rebound.Particle(m=1.)
    sim.add(sun)
    if wb_companion and companion_first:
        sim.add(primary=sun, m=1., a=500., e=0.2, inc=0.5, Omega=2., name="widebinary")
    sim.add(primary=sun, m=9.55e-4, a=5.2, e=0.05, inc=0.02, name="jupiter")
    sim.add(primary=sun, m=2.857e-4, a=9.58, e=0.06, inc=0.04, Omega=1., name="saturn")
    if wb_companion and not companion_first:
        sim.add(primary=sun, m=1., a=500., e=0.2, inc=0.5, Omega=2., name="widebinary")

    sim.move_to_com()
    sim.integrator = "trace"
    sim.integrator.coordinates = "widebinary"
    sim.dt = 2.*math.pi*5.2**1.5/20. # 1/20th Jupiter's orbital period
    return sim

def kozai_sim():
    sim = rebound.Simulation()
    sun = rebound.Particle(m=1.)
    sim.add(sun)
    sim.add(primary=sun,m=9.55e-4,a=5.2)
    sim.add(primary=sun,m=1,a=200,inc=math.pi/180.*75.,name="widebinary")
    return sim

class TestIntegratorTrace(unittest.TestCase):

    # ------ Boilerplate transformations/ordering ------------
    def test_transformation(self):
        clib = rebound.clibrebound
        # tests if transforming inertial -> WB -> inertial returns same ICs
        sim = jupsat_sim()
        before = sim.copy()
        clib.reb_integrator_trace_inertial_to_wb(ctypes.byref(sim))
        clib.reb_integrator_trace_wb_to_inertial(ctypes.byref(sim))
        after = sim.copy()

        for i in range(sim.N):
            self.assertAlmostEqual(before.particles[i].x, after.particles[i].x, delta=1e-8)
            self.assertAlmostEqual(before.particles[i].vx, after.particles[i].vx, delta=1e-8)
            self.assertAlmostEqual(before.particles[i].y, after.particles[i].y, delta=1e-8)
            self.assertAlmostEqual(before.particles[i].vy, after.particles[i].vy, delta=1e-8)
            self.assertAlmostEqual(before.particles[i].z, after.particles[i].z, delta=1e-8)
            self.assertAlmostEqual(before.particles[i].vz, after.particles[i].vz, delta=1e-8)

    def test_no_companion_matches_dhc(self):
        # Without a particle named "widebinary", WB coordinates should reduce to DHC.
        sim1 = jupsat_sim(wb_companion=False)
        sim2 = jupsat_sim(wb_companion=False)
        sim2.integrator.coordinates = "democraticheliocentric"
        with self.assertWarns(RuntimeWarning):
            sim1.steps(2000)
        sim2.steps(2000)
        for name in ["jupiter", "saturn"]:
            self.assertAlmostEqual(sim1.particles[name].x, sim2.particles[name].x, delta=1e-8)
            self.assertAlmostEqual(sim1.particles[name].vx, sim2.particles[name].vx, delta=1e-8)

    def test_companion_order_independence(self):
        # adding companion first or last should not change results
        sim1 = jupsat_sim(companion_first=False)
        sim2 = jupsat_sim(companion_first=True)
        sim1.steps(2000)
        sim2.steps(2000)
        for name in ["jupiter", "saturn", "widebinary"]:
            self.assertAlmostEqual(sim1.particles[name].x, sim2.particles[name].x, delta=1e-8)
            self.assertAlmostEqual(sim1.particles[name].vx, sim2.particles[name].vx, delta=1e-8)

    # ----- Test switching ----------------------------------------
    def test_no_switching_without_encounters(self):
        # checks for unnecessary pericenter switching when no encounters occur
        sim1 = jupsat_sim()
        sim2 = peri_switching_off(jupsat_sim())
        sim1.steps(2000)
        sim2.steps(2000)
        for name in ["jupiter", "saturn", "widebinary"]:
            self.assertAlmostEqual(sim1.particles[name].x, sim2.particles[name].x, delta=1e-13)
            self.assertAlmostEqual(sim1.particles[name].vx, sim2.particles[name].vx, delta=1e-13)

    def test_encounter_with_companion(self):
        # checks planet encounters with the wide binary companion
        def s_type_sim(integrator="trace"):
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(primary=sim.particles[0], m=0.5, a=30., name="widebinary")
            sim.add(primary=sim.particles["widebinary"], m=3e-6, a=0.5, name="planet") # within r_crit_WB of companion's Hill radius
            sim.move_to_com()
            sim.integrator = integrator
            if integrator == "trace":
                sim.integrator.coordinates = "widebinary"
            sim.dt = 0.1
            return sim
        
        switching = s_type_sim()
        no_switching = s_type_sim()
        no_switching.integrator.r_crit_WB = 0. # disable encounters with companion
        E0 = switching.energy()
        switching.integrate(30.)
        no_switching.integrate(30.)
        dE_switching = abs((switching.energy() - E0)/E0)
        dE_no_switching = abs((no_switching.energy() - E0)/E0)
    
        self.assertLess(dE_switching, 1e-11)                    # switching integrates the encounter accurately
        self.assertGreater(dE_no_switching, 1e3*dE_switching)   # and is far better than not switching

    
    def test_many_encounters(self):
        # tests planet-planet encounters
        def many_encounters_sim():
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
            sim.add(m=9.55e-4, a=5.2, name="widebinary")
            sim.move_to_com()
            sim.dt = 0.034
            return sim
        
        sim = many_encounters_sim()
        E0 = sim.energy()
        sim.integrate(500.)
        dE = abs((sim.energy() - E0)/E0)

        reference = many_encounters_sim()
        reference.integrator = "ias15"
        reference.integrate(500.)

        self.assertLess(dE, 1e-6)

    def test_no_effect_tp(self):
        # tests if test particle encounters have an effect
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=0.8,a=100.,name="widebinary")
        sim.move_to_com()
        sim.N_active = 3
        sim.integrator = "trace"
        sim.integrator.coordinates = 'widebinary'
        sim.dt = 0.1
        sim2 = sim.copy()
        p = sim.particles[1].copy()
        p.x += 0.01
        p.m = 0
        sim.add(p)
        sim.steps(1)
        sim2.steps(1)

        self.assertEqual(sim.particles[1].x,sim2.particles[1].x)
        self.assertEqual(sim.particles[1].vx,sim2.particles[1].vx)
        self.assertEqual(sim.particles[0].x,sim2.particles[0].x)
        self.assertEqual(sim.particles[0].vx,sim2.particles[0].vx)

    # ------ Restarts & additional forces

    def test_trace_simulationarchive(self):
        # Restarting from a snapshot must reproduce an uninterrupted run bit for bit.
        sim = chaotic_exchange_sim()
        sim.save_to_file("test.bin", step=10, delete_file=True)
        sim.integrate(1000., exact_finish_time=0)

        sim = rebound.Simulation("test.bin")
        sim.integrate(2000., exact_finish_time=0)
        x1 = sim.particles[1].x

        sim = chaotic_exchange_sim()
        sim.integrate(2000., exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0, x1)


    def test_additional_forces(self):
        # additional forces should be called twice per step
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.add(m=0.5, a=100., name="widebinary")
        sim.move_to_com()
        sim.integrator = "trace"
        sim.integrator.coordinates = "widebinary"
        sim.dt = 0.05

        calls = [0]
        def count_calls(reb_sim):
            calls[0] += 1
        sim.additional_forces = count_calls

        sim.steps(100)
        self.assertEqual(calls[0], 200)

    # ------ Collisions --------------

    def test_collision_with_star(self):
        sim = rebound.Simulation()
        sim.add(m=1.,r=0.00465)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1,f=2.3)
        sim.add(m=1e-4,r=1.4e-3,x=1.,vx=-0.4) # falling onto the star
        sim.add(m=1e-5,r=1.6e-4,a=1.5,e=0.1)
        sim.add(m=0.5,a=100.,name="widebinary")
        mtot0 = sum([p.m for p in sim.particles])
        com0 = sim.com()
        N0 = sim.N

        sim.integrator = "trace"
        sim.integrator.coordinates = 'widebinary'
        sim.dt = 0.01
        sim.track_energy_offset = 1
        sim.collision = "direct"
        sim.collision_resolve = "merge"

        E0 = sim.energy()
        sim.integrate(1)
        com1 = sim.com()
        self.assertAlmostEqual(com1.vx,com0.vx,delta=2e-16)
        self.assertAlmostEqual(com1.vy,com0.vy,delta=2e-16)
        self.assertAlmostEqual(com1.vz,com0.vz,delta=2e-16)
        mtot1 = sum([p.m for p in sim.particles])
        self.assertEqual(mtot0,mtot1)
        self.assertEqual(N0-1,sim.N)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,5e-2)

    def test_collision_add_particles(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-5,r=1.6e-4,a=0.5,e=0.1)    #these params lead to collision on my machine
        sim.add(m=1e-8,r=4e-5,a=0.55,e=0.4,f=-0.94)
        sim.add(m=0.5,a=50.,name="widebinary")
        N0 = sim.N

        sim.integrator = "trace"
        sim.integrator.coordinates = 'widebinary'
        sim.dt = 0.01
        sim.collision = "direct"
        sim.collision_resolve = collision_add_particle

        sim.integrate(1)
        self.assertEqual(N0+1,sim.N)

        E0 = sim.energy()
        sim.integrate(1000.)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,1e-3)

    # --------- Accuracy Comparisons -----------------

    def outer_solar_system(self):
        def max_energy_error(integrator):
            sim = rebound.Simulation()
            sim.add("outer solar system")
            sim.add(m=0.5, a=200., name="widebinary")
            sim.integrator = integrator
            if integrator == "trace":
                sim.integrator.coordinates = "widebinary"
            sim.dt = sim.particles[1].P/20.
            E0 = sim.energy()

            # calculate maximum relative energy error
            dE_max = 0.
            for k in range(1, 201):
                sim.integrate(5000.*k/200., exact_finish_time=0)
                dE_max = max(dE_max, abs((sim.energy() - E0)/E0))
            return dE_max
        dE_wb = max_energy_error("trace")
        dE_whfast = max_energy_error("whfast")
        self.assertLess(dE_wb, 5.*dE_whfast) # competitive with WHFast

    def test_kozai(self):
        def max_eccentricity(integrator):
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(primary=sim.particles[0], m=9.55e-4, a=5.2, e=0.05, name="jupiter")
            sim.add(primary=sim.particles[0], m=1., a=50., inc=75.*math.pi/180., name="widebinary")
            sim.move_to_com()
            sim.integrator = integrator
            if integrator == "trace":
                sim.integrator.coordinates = "widebinary"
            sim.dt = 2.*math.pi*5.2**1.5/30.
            t_zlk = 16./(30.*math.pi)*(2. + 9.55e-4)*(2.*math.pi*math.sqrt(50.**3/2.))**2/(2.*math.pi*5.2**1.5)
            e_max = 0.
            for k in range(1, 301):
                sim.integrate(4.*t_zlk*k/300., exact_finish_time=0)
                e_max = max(e_max, sim.particles["jupiter"].orbit(primary=sim.particles[0]).e)
            return e_max
        e_max_wb = max_eccentricity("trace")
        e_max_ias15 = max_eccentricity("ias15")
        self.assertGreater(e_max_ias15, 0.9)                  # the ZLK cycle reaches high eccentricity
        self.assertLess(abs(e_max_wb - e_max_ias15), 5e-3)

    def test_outer_solar_massive(self):
        sim = rebound.Simulation()
        sim.add("outer solar system")
        for i in range(1,sim.N):
            sim.particles[i].m *=50.
        sim.add(m=1.,a=200., name="widebinary")

        sim.integrator = "trace"
        sim.integrator.coordinates = 'widebinary'
        P = sim.particles[1].P
        sim.dt = 1e-3*P

        E0 = sim.energy()
        sim.integrate(1000)
        dE = abs((sim.energy() - E0)/E0)
        self.assertLess(dE,1e-6)

if __name__ == "__main__":
    unittest.main()
