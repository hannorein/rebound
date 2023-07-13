import rebound
import unittest
import os
import numpy as np
import time

def GetApophis1979():
    return (np.array([[ 2.0461278699861521e-06,  2.2370897738299570e-06,
             -1.6661854822847768e-10],
            [-6.8124977762585481e-01, -7.4482975052825995e-01,
              5.5474953737657817e-05],
            [ 7.5476934423654385e-01, -1.8985795987726539e-03,
              1.8667813859384475e-02]]),
     np.array([[-3.7290167155085508e-08,  3.5092072454320220e-08,
             -2.6790975184377306e-12],
            [ 1.2415606304314073e-02, -1.1683759801433082e-02,
              8.9199439363293445e-07],
            [ 1.7223306120259661e-03,  2.1389341162818411e-02,
             -1.0892481215440922e-03]]),
     np.array([2.959122082855911e-04, 8.887697821383033e-10,
            4.016599676151095e-24]),
     365.24141698304345,
     False)    

def create_sim():
    sim = rebound.Simulation()
    sim.G = 1.4881806877180788e-34
    problem = GetApophis1979
    Q,V,mass,period,_ = problem()
    mass /= sim.G
    
    for i in range(3):
        sim.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
    
    sim.move_to_com()
    sim.dt = period/100
    sim.integrator = "tes"  
    return sim, period

class TestIntegratorTES(unittest.TestCase):
    def test_integration_output_particles(self):  
        orbits = 100
        problem = GetApophis1979
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
        sim.ri_tes.orbits = orbits
                
        e0 = sim.energy()
        sim.integrate(period*orbits)
        e1 = sim.energy()
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
        problem = GetApophis1979
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
        problem = GetApophis1979
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
                
        e0 = sim.energy()
        times = np.linspace(0, orbits*period, 10)
        for t in times:
            sim.integrate(t)
        e1 = sim.energy()
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
        self.assertLess(de, 2e-15)
        
    def test_add_remove_particle(self):  
        orbits = 1
        problem = GetApophis1979
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
                
        e0 = sim.energy()
        outputs=100
        times = np.linspace(0, orbits*period, outputs)
        pos_out = np.zeros([outputs, 2, 3])
        for i, t in enumerate(times):
            sim.integrate(t)
            for j in range(2):
                pos_out[i,j,0] = sim.particles[j].x
                pos_out[i,j,1] = sim.particles[j].y
                pos_out[i,j,2] = sim.particles[j].z
                
        e1 = sim.energy()
        de = np.abs((e1-e0)/e0)
        self.assertLess(de, 1e-15)
        
        
        # Add Apophis to the simulation.
        sim.add(m=mass[2], x=Q[2,0], y=Q[2,1], z=Q[2,2], vx=V[2,0], vy=V[2,1], vz=V[2,2], hash=2)
        e0 = sim.energy()
        outputs=100
        times = np.linspace(orbits*period, 2*orbits*period, outputs)
        pos_out = np.zeros([outputs, 3, 3])
        for i, t in enumerate(times):
            sim.integrate(t)
            for j in range(3):
                pos_out[i,j,0] = sim.particles[j].x
                pos_out[i,j,1] = sim.particles[j].y
                pos_out[i,j,2] = sim.particles[j].z
                
        e1 = sim.energy()
        de = np.abs((e1-e0)/e0)        
        self.assertLess(de, 1e-15)
        
        
        # Remove Apophis from the simulation.
        sim.remove(2)
        e0 = sim.energy()
        outputs=100
        times = np.linspace(2*orbits*period, 3*orbits*period, outputs)
        pos_out = np.zeros([outputs, 3, 3])
        for i, t in enumerate(times):
            sim.integrate(t)
            for j in range(2):
                pos_out[i,j,0] = sim.particles[j].x
                pos_out[i,j,1] = sim.particles[j].y
                pos_out[i,j,2] = sim.particles[j].z
                
        e1 = sim.energy()
        de = np.abs((e1-e0)/e0)        
        self.assertLess(de, 1e-15)        
        
        
    def test_simulation_archive_bitwise(self):                        
        # sim0 is the ground truth
        sim0, period = create_sim()
        sim0.integrate(100*period, exact_finish_time=1)
        
        try:
            os.remove("test_simulation_archive_bitwise.archive")        
        except:
            pass
        
        # Perform half an integration and save archive
        sim1,_ = create_sim()
        sim1.integrate(50*period, exact_finish_time=0)
        sim1.simulationarchive_snapshot("test_simulation_archive_bitwise.archive")
        
        # Load archive and do the second half of the integration
        sim2 = rebound.SimulationArchive("test_simulation_archive_bitwise.archive")[-1]
        sim2.integrate(100*period, exact_finish_time=1)
        
        self.assertEqual(sim0.particles[1].x, sim2.particles[1].x)
        self.assertEqual(sim0.particles[1].y, sim2.particles[1].y)
        self.assertEqual(sim0.particles[1].z, sim2.particles[1].z)
        
        try:
            os.remove("test_simulation_archive_bitwise.archive")        
        except:
            pass
        
    def test_add_particle_then_save(self):                        
        # sim0 is the ground truth
        sim0, period = create_sim()
        sim0.integrate(100*period, exact_finish_time=1)
        
        try:
            os.remove("test_simulation_archive_bitwise.archive")        
        except:
            pass
        
        # Perform half an integration and save archive
        sim1,_ = create_sim()
        sim1.integrate(50*period, exact_finish_time=0)
        sim1.add(m=1, x=1.0, y=1.1, z=1.2, vx=0.0, vy=1.0, vz=0.0)
        sim1.simulationarchive_snapshot("test_simulation_archive_bitwise.archive")

        # Load archive and check for extra particle
        sim2 = rebound.SimulationArchive("test_simulation_archive_bitwise.archive")[-1]
        
        self.assertEqual(sim2.particles[-1].x, 1.0)
        self.assertEqual(sim2.particles[-1].y, 1.1)
        self.assertEqual(sim2.particles[-1].z, 1.2)
        
        try:
            os.remove("test_simulation_archive_bitwise.archive")        
        except:
            pass        
            
    def test_longterm_conservation(self): 
        return # Comment to run test (takes a few minutes)
        Q = np.array([[ 1.6166089996756736e-07, -1.3174811284580585e-06,  2.7324913024933070e-08],
                       [ 2.5649647880625376e-01,  1.9231594858567522e-01, -7.8597568865791164e-03],
                       [-3.5944957551650537e-02, -7.2591833677411899e-01, -7.7858351034036984e-03],
                       [-1.8053334772678853e-01,  9.6659467393463594e-01, 6.4783933279786331e-05],
                       [ 1.3261049612058200e+00,  4.9619050870453013e-01, -2.2292003930474757e-02]])

        V = np.array([[ 7.6312215682912423e-09,  3.7792298576194310e-09,  2.0802159000845027e-09],
                   [-2.2401743922574373e-02,  2.3733227399426897e-02,     3.9957465111722101e-03],
                   [ 2.0064892939487048e-02, -1.0775553592363494e-03,    -1.1732906261155678e-03],
                   [-1.7193861432114162e-02, -3.2158096004801006e-03,    2.3079228394344072e-07],
                   [-4.3648866780331211e-03,  1.4300831106604298e-02,   4.0702686654480900e-04]])

        mass = np.array([2.959122082855911e-04, 4.888673559153889e-11, 7.242975407123889e-10, 8.887415067052367e-10, 9.509474594518523e-11])
        period = 87.96828365969395    
        
        sim = rebound.Simulation()
        
        sim.G = 1.4881806877180788e-34
        mass /= sim.G
        
        for i in range(3):
            sim.add(m=mass[i], x=Q[i,0], y=Q[i,1], z=Q[i,2], vx=V[i,0], vy=V[i,1], vz=V[i,2])
        
        sim.move_to_com()
        # sim.dt = period/100
        sim.integrator = "tes"          
        e0 = sim.energy()

        samples = 10000
        orbits = 100
        e_arr = np.zeros(samples)
        t_arr = np.zeros(samples)
        for i in range(samples):
            sim.integrate(i*orbits*period)
            e_arr[i] = sim.energy()
            t_arr[i] = i*orbits
        de = np.abs((e_arr-e0)/e0)
        from matplotlib import pyplot as plt
        plt.figure(dpi=300)
        plt.loglog(t_arr[1:], de[1:])
        plt.ylim(1e-16, 1e-13)
        plt.grid()
        plt.xlabel("time (orbits)")
        plt.ylabel("relative energy error")
    
if __name__ == "__main__":
    unittest.main()

    

    
