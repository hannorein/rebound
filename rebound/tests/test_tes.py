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

class TestIntegratorTES(unittest.TestCase):
    def test_particle_pass_through(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3,a=1.12313)
        sim.add(m=1e-3,a=2.32323)
        sim.move_to_com()
        sim.dt = 0.25
        sim.integrator = "tes"
        pos = sim.particles[1].x
        sim.integrate(1)

        self.assertEqual(pos, sim.particles[1].x)

    # Generic test case for energy - yet to be implemented. 
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
        
        e0 = sim.calculate_energy()
        sim.integrate(1)
        e1 = sim.calculate_energy() 
    
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
        
        
if __name__ == "__main__":
    unittest.main()
    

    
    
    

    