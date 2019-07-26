import rebound
import unittest
import math
import rebound.data
    
sabasettings1 = [ # k, corrector, relative error
        [2,"none",1e-10], 
        [2,"modifiedkick",1e-10], 
        [2,"lazy",1e-10], 
        [3,"none",2e-11], 
        [3,"modifiedkick",2e-11], 
        [3,"lazy",2e-11], 
        [4,"none",2e-11], 
        [4,"modifiedkick",2e-11], 
        [4,"lazy",2e-11], 
        ]
sabasettings2 = [ # k, corrector, relative error
        [1,"none",2e-14], 
        [1,"modifiedkick",4e-14], 
        [1,"lazy",4e-14], 
        [2,"none",4e-14], 
        [2,"modifiedkick",4e-14], 
        [2,"lazy",4e-14], 
        [3,"none",4e-14], 
        [3,"modifiedkick",4e-14], 
        [3,"lazy",4e-14], 
        [4,"none",4e-14], 
        [4,"modifiedkick",4e-14], 
        [4,"lazy",4e-14], 
        ]

class TestIntegratorSABA(unittest.TestCase):
    def test_sabasettings1(self):
        for s in sabasettings1:
            test_name = "test_energy_SABA_k_%02d_c_%s" % (s[0], s[1])
            self.energy(s)
            test_name = "test_energy_notcom_SABA_k_%02d_c_%s" % (s[0], s[1])
            self.energy_notcom(s)
            test_name = "test_compias_SABA_k_%02d_c_%s" % (s[0], s[1])
            self.compias(s)
    def test_sabasettings2(self):
        for s in sabasettings2:
            test_name = "test_backandforth_SABA_k_%02d_c_%s" % (s[0], s[1])
            self.backandforth(s)
    def test_sabarestart(self):
        for s in sabasettings2:
            test_name = "test_restart_SABA_k_%02d_c_%s" % (s[0], s[1])
            self.restart(s)

    def energy(self, s):
        k, corrector, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.integrator = "saba"
        sim.ri_saba.corrector = corrector 
        sim.ri_saba.k = k
        sim.ri_saba.safe_mode = False
        sim.dt = 0.0123235235*sim.particles[1].P  
        e0 = sim.calculate_energy()
        sim.integrate(1000.*2.*3.1415)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),maxerror)

    def energy_notcom(self, s):
        k, corrector, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for p in sim.particles:
            p.vx += 1.
        com = sim.calculate_com()
        sim.integrator = "saba"
        sim.ri_saba.corrector = corrector 
        sim.ri_saba.k = k
        sim.dt = 0.0123235235*sim.particles[1].P  
        e0 = sim.calculate_energy()
        sim.integrate(1000.*2.*3.1415)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),maxerror)
        com1 = sim.calculate_com()
        self.assertLess(math.fabs((com.x+com.vx*sim.t-com1.x)/(com1.x+com1.y)),1e-12)
        self.assertLess(math.fabs((com.y+com.vy*sim.t-com1.y)/(com1.x+com1.y)),1e-12)

    def compias(self, s):
        k, corrector, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for p in sim.particles:
            p.vx += 1050. # move out of com to make it harder
        sim.integrator = "saba"
        sim.ri_saba.corrector = corrector 
        sim.ri_saba.k = k
        sim.dt = 0.0123235235*sim.particles[1].P  
        sim.integrate(13.21415,exact_finish_time=False)
        
        simi = rebound.Simulation()
        simi.integrator = "ias15"
        rebound.data.add_outer_solar_system(simi)
        for p in simi.particles:
            p.vx += 1050.
        simi.integrate(sim.t,exact_finish_time=True)
        
        for i in range(sim.N):
            self.assertLess(math.fabs(simi.particles[i].x-sim.particles[i].x),5e-9)

    def backandforth(self, s):
        k, corrector, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for p in sim.particles:
            p.vx += 1.
        sim0=sim.copy()
        sim.integrator = "saba"
        sim.ri_saba.corrector = corrector 
        sim.ri_saba.k = k
        sim.dt = 0.0123235235*sim.particles[1].P  
        steps = 10
        for i in range(steps):
            sim.step()
        sim.dt *= -1
        for i in range(steps):
            sim.step()
        for i in range(sim.N):
            self.assertLess(math.fabs(sim0.particles[i].x-sim.particles[i].x),maxerror)
    
    def restart(self, s):
        k, corrector, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.integrator = "saba"
        sim.ri_saba.corrector = corrector 
        sim.ri_saba.k = k
        sim.step()
        sim2 = sim.copy()
        sim.step()
        sim2.step()
        self.assertEqual(sim,sim2)

if __name__ == "__main__":
    unittest.main()
