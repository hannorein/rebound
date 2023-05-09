import rebound
import unittest
import math
import rebound.data
import warnings

class TestIAS15mode1(unittest.TestCase):
    
    def test_safe_and_load(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.move_to_com()
        sim.integrator = "ias15"
        sim.ri_ias15.dt_mode = 1
        sim.save("sim.bin")
        jupyr = 11.86*2.*math.pi
        sim.integrate(1e3*jupyr)

        sim2 = rebound.Simulation("sim.bin")
        sim2.integrate(-1e3*jupyr)
        
        for i in range(sim.N):
            self.assertEqual(sim2.particles[i].x,sim.particles[i].x)
            self.assertEqual(sim2.particles[i].y,sim.particles[i].y)
            self.assertEqual(sim2.particles[i].z,sim.particles[i].z)
    
    def test_backwards(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.move_to_com()
        sim.integrator = "ias15"
        sim.ri_ias15.dt_mode = 1
        sim2 = sim.copy()
        jupyr = 11.86*2.*math.pi
        sim.integrate(1e3*jupyr)

        for p in sim2.particles:
            p.vx *= -1
            p.vy *= -1
            p.vz *= -1
        sim2.dt *= -1
        sim2.integrate(-1e3*jupyr)
        
        for i in range(sim.N):
            self.assertEqual(sim2.particles[i].x,sim.particles[i].x)
            self.assertEqual(sim2.particles[i].y,sim.particles[i].y)
            self.assertEqual(sim2.particles[i].z,sim.particles[i].z)
    
    def test_outer_solar_system(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        sim.move_to_com()
        sim.integrator = "ias15"
        sim.ri_ias15.dt_mode = 1
        jupyr = 11.86*2.*math.pi
        e0 = sim.energy()
        self.assertNotEqual(e0,0.)
        sim.integrate(1e3*jupyr)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e1),1e-15)


if __name__ == "__main__":
    unittest.main()
