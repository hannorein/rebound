import rebound
import unittest
import math
import rebound.data
    
sabasettings1 = [ # type, relative error
        ["SABA2",1e-10], 
        ["SABACM2",1e-10], 
        ["SABACL2",1e-10], 
        ["SABA3",2e-11], 
        ["SABACM3",2e-11], 
        ["SABACL3",2e-11], 
        ["SABA4",2e-11], 
        ["SABACM4",2e-11], 
        ["SABACL4",2e-11], 
        ["SABA10,4",5e-13], 
        ["SABA8,6,4",5e-13], 
        ["SABA10,6,4",5e-13], 
        ["SABAH8,4,4",5e-13], 
        ["SABAH8,6,4",5e-13], 
        ["SABAH10,6,4",5e-13], 
        ]
sabasettings2 = [ # type, relative error
        ["SABA1",2e-14], 
        ["SABACM4",4e-14], 
        ["SABACL4",4e-14], 
        ["SABA2",4e-14], 
        ["SABACM4",4e-14], 
        ["SABACL4",4e-14], 
        ["SABA3",4e-14], 
        ["SABACM4",4e-14], 
        ["SABACL4",4e-14], 
        ["SABA4",4e-14], 
        ["SABACM4",4e-14], 
        ["SABACL4",4e-14], 
        ["SABA10,4",3e-14], 
        ["SABA8,6,4",3e-14], 
        ["SABA10,6,4",3e-14], 
        ["SABAH8,4,4",3e-14], 
        ["SABAH8,6,4",3e-14], 
        ["SABAH10,6,4",3e-14], 
        ]

class TestIntegratorSABA(unittest.TestCase):
    def test_sabasettings1(self):
        for s in sabasettings1:
            test_name = "test_energy_%s" % (s[0])
            self.energy(s)
            test_name = "test_energy_notcom_%s" % (s[0])
            self.energy_notcom(s)
            test_name = "test_compias_%s" % (s[0])
            self.compias(s)
    def test_sabasettings2(self):
        for s in sabasettings2:
            test_name = "test_backandforth_%s" % (s[0])
            self.backandforth(s)
    def test_sabarestart(self):
        for s in sabasettings2:
            test_name = "test_restart_%s" % (s[0])
            self.restart(s)

    def energy(self, s):
        integrator, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.integrator = integrator
        sim.ri_saba.safe_mode = False
        sim.dt = 0.0123235235*sim.particles[1].P  
        e0 = sim.calculate_energy()
        sim.integrate(1000.*2.*3.1415)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),maxerror)

    def energy_notcom(self, s):
        integrator, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for p in sim.particles:
            p.vx += 1.
        com = sim.calculate_com()
        sim.integrator = integrator
        sim.dt = 0.0123235235*sim.particles[1].P  
        e0 = sim.calculate_energy()
        sim.integrate(1000.*2.*3.1415)
        e1 = sim.calculate_energy()
        self.assertLess(math.fabs((e0-e1)/e1),maxerror)
        com1 = sim.calculate_com()
        self.assertLess(math.fabs((com.x+com.vx*sim.t-com1.x)/(com1.x+com1.y)),1e-12)
        self.assertLess(math.fabs((com.y+com.vy*sim.t-com1.y)/(com1.x+com1.y)),1e-12)

    def compias(self, s):
        integrator, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for p in sim.particles:
            p.vx += 1050. # move out of com to make it harder
        sim.integrator = integrator
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
        integrator, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        for p in sim.particles:
            p.vx += 1.
        sim0=sim.copy()
        sim.integrator = integrator
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
        integrator, maxerror = s
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.integrator = integrator
        sim.step()
        sim2 = sim.copy()
        sim.step()
        sim2.step()
        self.assertEqual(sim,sim2)

if __name__ == "__main__":
    unittest.main()
