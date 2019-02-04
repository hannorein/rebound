import rebound
import unittest
import math
import rebound.data
import warnings

class TestIntegratorJanus(unittest.TestCase):
    def test_janus_energy(self):
        for o, eps in [ (2,1e-4), (4,1e-8), (6,1e-9), (8,1e-11), (10,1e-13)]:
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(m=1e-3,a=1.12313)
            sim.add(m=1e-3,a=2.32323)
            sim.move_to_com()
            sim.dt = 0.25
            sim.integrator = "janus"
            sim.ri_janus.order = o
            sim.ri_janus.scale_pos = 1e-16
            sim.ri_janus.scale_vel = 1e-16
            e0 = sim.calculate_energy()
            sim.integrate(1e2)
            e1 = sim.calculate_energy()
            self.assertLess(math.fabs((e0-e1)/e1),eps)
    
    def test_janus_reverse(self):
        for o in [2,4,6,8,10]:
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(m=1e-3,a=1.12313,omega=0.32643,l=0.3788,e=0.012)
            sim.add(m=1e-3,a=2.32323,omega=0.12314,l=0.1726,e=0.103)
            sim.move_to_com()
            sim.dt = 0.25
            sim.integrator = "janus"
            sim.ri_janus.order = o
            sim.ri_janus.safe_mode = 0
            sim.ri_janus.scale_pos = 1e-16
            sim.ri_janus.scale_vel = 1e-16
            sim.step()
            t0 = sim.t

            x1, x2 = sim.particles[1].x, sim.particles[2].x
            vx1, vx2 = sim.particles[1].vx, sim.particles[2].vx
            y1, y2 = sim.particles[1].y, sim.particles[2].y
            vy1, vy2 = sim.particles[1].vy, sim.particles[2].vy

            sim.integrate(1e2,exact_finish_time=0)
            sim.dt *= -1
            sim.integrate(t0,exact_finish_time=0)
            
            xf1, xf2 = sim.particles[1].x, sim.particles[2].x
            vxf1, vxf2 = sim.particles[1].vx, sim.particles[2].vx
            yf1, yf2 = sim.particles[1].y, sim.particles[2].y
            vyf1, vyf2 = sim.particles[1].vy, sim.particles[2].vy
            
            self.assertEqual(x1,xf1)
            self.assertEqual(x2,xf2)
            self.assertEqual(vx1,vxf1)
            self.assertEqual(vx2,vxf2)
            self.assertEqual(y1,yf1)
            self.assertEqual(y2,yf2)
            self.assertEqual(vy1,vyf1)
            self.assertEqual(vy2,vyf2)
    
    def test_janus_simulationarchive(self):
        for o in [2,4,6,8,10]:
            sim = rebound.Simulation()
            sim.add(m=1.)
            sim.add(m=1e-3,a=1.12313)
            sim.add(m=1e-3,a=2.32323)
            sim.move_to_com()
            sim.dt = 0.25
            sim.automateSimulationArchive("test.bin",interval=5,deletefile=True)
            sim.integrator = "janus"
            sim.ri_janus.order = o
            sim.ri_janus.scale_pos = 1e-16
            sim.ri_janus.scale_vel = 1e-16
            sim.integrate(1e2,exact_finish_time=0)
            
            sim2 = rebound.Simulation("test.bin")
            sim2.integrate(2e2,exact_finish_time=0)
            
            sim.integrate(2e2,exact_finish_time=0)

            self.assertEqual(sim.t,sim2.t);
            self.assertEqual(sim.particles[1].x,sim2.particles[1].x);
            
            


if __name__ == "__main__":
    unittest.main()
