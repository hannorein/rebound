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
    # Test case for energy - yet to be implemented. 
    # def test_particle_pass_through(self):
        # sim = rebound.Simulation()
        # sim.add(m=1.)
        # sim.add(m=1e-3,a=1.12313)
        # sim.add(m=1e-3,a=2.32323)
        # sim.move_to_com()
        # sim.dt = 0.5
        # sim.integrator = "tes"
        # e0 = sim.calculate_energy()
        # sim.integrate(1)
        # e1 = sim.calculate_energy()
        # self.assertLess(math.fabs((e0-e1)/e1),eps)


    def test_bitwise_identical_tes_v1_vs_tes_rebound(self):
        orbits = 100
        problem = init.GetApophis1979
        output_samples=2500
        samples = 1
        tol=1e-6
        recti_per_orbit = 1.61803398875    
        
        ###################################################
        # Run experiments using TES in REBOUND
        ###################################################       
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
        
        sim.integrate(1)    
        data = tes_driver.extract_data('temp_output.txt')
    
    
        ###################################################
        # Run experiments using TES v1 distibution
        ###################################################
    
        manager = Exp.ExperimentManager("apophis_baseline") 
        config = {}
        config['orbits'] = orbits
        config['rtol'] = tol
        config['name'] = 'tes'
        config['integrator'] = 'tes'
        config['output samples'] = output_samples
        config['output spacing'] = 'linear'   
        config['dQ max'] = 1E-3
        config['dP max'] = 1
        config['rectifications per orbit'] = 1.61803398875
        config['steps initial per orbit'] = 1E2
        config['time out'] = 2*86400
        config['exe file'] = 'output_v1'
        manager.add_experiment(config, problem)
       
        manager.run_experiments()
        manager.save_experiments()    
      
        
        ###################################################
        # Comparison between the two tools
        ###################################################
        data_reb = manager.experiments[0].integration
        data_tes = data
        
        error = np.abs(data_tes['Q'] - data_reb['Q'])
        error = error[error > 0.0]
        self.assertEqual(len(error), 0)
        
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
        
        sim.integrate(period*orbits)
        
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
        
        error = np.abs(data_tes - tes_reb_pos)
        errors = len(error[error > 0.0])
        self.assertEqual(errors, 0)                    
        
    def test_timing_with_ias15(self):       
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
        self.assertGreater(rt_ias, rt_tes)  


if __name__ == "__main__":
    unittest.main()

    
    

    
