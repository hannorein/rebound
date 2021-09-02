import rebound
import unittest
import math
import rebound.data
import warnings
import copy
import os
import sys
import pathlib
sys.path.append('/home/pete/Desktop/TES_V1/experiment_manager/experiments/')
sys.path.append('/home/pete/Desktop/TES_V1/experiment_manager/drivers/')
from matplotlib import pyplot as plt
import initial_conditions as init
import experiment_manager as Exp
import numpy as np
import tes_driver
import time

class TestIntegratorTES(unittest.TestCase):
    def test_integration_output_particles(self):  
        orbits = 100
        problem = init.GetApophis1979
        output_samples=2500
        samples = 1
        tol=1e-6
        recti_per_orbit = 1.61803398875    
                 
        G_au_kg_dy = 1.4881806877180788e-34   
        Q,V,mass,period,_ = problem()
        mass /= G_au_kg_dy
        
        sim = rebound.Simulation()
        sim.G = G_au_kg_dy
        for i in range(3):
            sim.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
        
        sim.move_to_com()
        sim.dt = period/100
        sim.integrator = "tes"  
        sim.ri_tes.dq_max = 1e-3
        sim.ri_tes.recti_per_orbit = recti_per_orbit
        sim.ri_tes.epsilon = tol
        sim.ri_tes.output_samples = output_samples
        sim.ri_tes.orbital_period = period
        sim.ri_tes.orbits = orbits
                
        e0 = sim.calculate_energy()
        sim.integrate(period*orbits)
        e1 = sim.calculate_energy()
        de = np.abs((e1-e0)/e0)

        data_tes = np.array([[ 0.0000000000000000e+00,  0.0000000000000000e+00,
             0.0000000000000000e+00],
           [-9.9604746697684021e-01,  3.5311915613020404e-03,
            -1.2054180564475472e-06],
           [-8.1089946547081804e-01, -5.4094730893500276e-01,
             6.8972157890442951e-03]])
    
        N=3
        tes_reb_pos = np.zeros([N,3])
        for i in range(N):
            tes_reb_pos[i,0] = sim.particles[i].x
            tes_reb_pos[i,1] = sim.particles[i].y
            tes_reb_pos[i,2] = sim.particles[i].z    
        
        error = np.max(np.abs(data_tes - tes_reb_pos))
        # 1e-7 is max precision (from article) - noise here is due to dh coords change.
        self.assertLess(error, 5e-6) 
        self.assertLess(de, 1e-15)
        
        
    def test_timing_with_ias15(self):    
        orbits = 100
        problem = init.GetApophis1979
        output_samples=1
        samples = 1
        tol=1e-6
        recti_per_orbit = 1.61803398875    
                 
        G_au_kg_dy = 1.4881806877180788e-34   
        Q,V,mass,period,_ = problem()
        mass /= G_au_kg_dy
        
        sim = rebound.Simulation()
        sim.G = G_au_kg_dy
        for i in range(3):
            sim.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
        
        sim.move_to_com()
        sim.dt = period/100
        sim.integrator = "tes"  
        sim.ri_tes.dq_max = 1e-3
        sim.ri_tes.recti_per_orbit = recti_per_orbit
        sim.ri_tes.epsilon = tol
        sim.ri_tes.output_samples = output_samples
        sim.ri_tes.orbital_period = period
        sim.ri_tes.orbits = orbits
        
        t0_tes = time.time()
        sim.integrate(period*orbits)
        t1_tes = time.time()
      
        sim2 = rebound.Simulation()
        sim2.G = G_au_kg_dy
        for i in range(3):
            sim2.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
        
        sim2.move_to_com()
        sim2.integrator = "ias15"  
        
        t0_ias = time.time()
        sim2.integrate(period*orbits)
        t1_ias = time.time()

        rt_ias = t1_ias-t0_ias
        rt_tes = t1_tes-t0_tes
        print('\nruntime tes/ias15:', rt_tes/rt_ias)
        self.assertGreater(rt_ias, rt_tes)  

    def test_piecewise_integration(self):  
        orbits = 100
        problem = init.GetApophis1979
        output_samples=2500
        samples = 1
        tol=1e-6
        recti_per_orbit = 1.61803398875  
                 
        G_au_kg_dy = 1.4881806877180788e-34   
        Q,V,mass,period,_ = problem()
        mass /= G_au_kg_dy
        
        sim = rebound.Simulation()
        sim.G = G_au_kg_dy
        for i in range(3):
            sim.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
        
        sim.move_to_com()
        sim.dt = period/100
        sim.integrator = "tes"  
        sim.ri_tes.dq_max = 1e-3
        sim.ri_tes.recti_per_orbit = recti_per_orbit
        sim.ri_tes.epsilon = tol
                
        e0 = sim.calculate_energy()
        times = np.linspace(0, orbits*period, 10)
        for t in times:
            sim.integrate(t)
        e1 = sim.calculate_energy()
        de = np.abs((e1-e0)/e0)

        data_tes = np.array([[ 0.0000000000000000e+00,  0.0000000000000000e+00,
             0.0000000000000000e+00],
           [-9.9604746697684021e-01,  3.5311915613020404e-03,
            -1.2054180564475472e-06],
           [-8.1089946547081804e-01, -5.4094730893500276e-01,
             6.8972157890442951e-03]])
    
        N=3
        tes_reb_pos = np.zeros([N,3])
        for i in range(N):
            tes_reb_pos[i,0] = sim.particles[i].x
            tes_reb_pos[i,1] = sim.particles[i].y
            tes_reb_pos[i,2] = sim.particles[i].z    
        
        error = np.max(np.abs(data_tes - tes_reb_pos))
        # 1e-7 is max precision (from article) - noise here is due to dh coords change.
        self.assertLess(error, 5e-6) 
        self.assertLess(de, 1e-15)
        
    def test_add_remove_particle(self):  
        orbits = 1
        problem = init.GetApophis1979
        output_samples=2500
        samples = 1
        tol=1e-6
        recti_per_orbit = 1.61803398875  
                 
        G_au_kg_dy = 1.4881806877180788e-34   
        Q,V,mass,period,_ = problem()
        mass /= G_au_kg_dy
        
        sim = rebound.Simulation()
        sim.G = G_au_kg_dy
        
        # Add the sun and Earth
        for i in range(2):
            sim.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
        
        sim.move_to_com()
        sim.dt = period/100
        sim.integrator = "tes"  
        sim.ri_tes.dq_max = 1e-3
        sim.ri_tes.recti_per_orbit = recti_per_orbit
        sim.ri_tes.epsilon = tol
                
        e0 = sim.calculate_energy()
        outputs=100
        times = np.linspace(0, orbits*period, outputs)
        pos_out = np.zeros([outputs, 2, 3])
        for i, t in enumerate(times):
            sim.integrate(t)
            for j in range(2):
                pos_out[i,j,0] = sim.particles[j].x
                pos_out[i,j,1] = sim.particles[j].y
                pos_out[i,j,2] = sim.particles[j].z
                
        e1 = sim.calculate_energy()
        de = np.abs((e1-e0)/e0)
        self.assertLess(de, 1e-15)
        
        
        # Add Apophis to the simulation.
        sim.add(m=mass[2], x=Q[2,0], y=Q[2,1], z=Q[2,2], vx=V[2,0], vy=V[2,1], vz=V[2,2], hash=2)
        e0 = sim.calculate_energy()
        outputs=100
        times = np.linspace(orbits*period, 2*orbits*period, outputs)
        pos_out = np.zeros([outputs, 3, 3])
        for i, t in enumerate(times):
            sim.integrate(t)
            for j in range(3):
                pos_out[i,j,0] = sim.particles[j].x
                pos_out[i,j,1] = sim.particles[j].y
                pos_out[i,j,2] = sim.particles[j].z
                
        e1 = sim.calculate_energy()
        de = np.abs((e1-e0)/e0)        
        self.assertLess(de, 1e-15)
        
        
        # Remove Apophis from the simulation.
        sim.remove(sim.particles[2])
        e0 = sim.calculate_energy()
        outputs=100
        times = np.linspace(2*orbits*period, 3*orbits*period, outputs)
        pos_out = np.zeros([outputs, 3, 3])
        for i, t in enumerate(times):
            sim.integrate(t)
            for j in range(2):
                pos_out[i,j,0] = sim.particles[j].x
                pos_out[i,j,1] = sim.particles[j].y
                pos_out[i,j,2] = sim.particles[j].z
                
        e1 = sim.calculate_energy()
        de = np.abs((e1-e0)/e0)        
        self.assertLess(de, 1e-15)        
        
if __name__ == "__main__":
    unittest.main()

    
    

    
