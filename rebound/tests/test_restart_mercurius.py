import rebound
import unittest
import warnings

class TestSimulationRestartMercurius(unittest.TestCase):
    def test_sa_mercurius_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "mercurius"
        sim.dt = 0.1313
        sim.ri_mercurius.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)
        sim.save("test.bin")
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        sim = rebound.Simulation.from_file("test.bin")
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    def test_sa_mercurius_restart_safemode(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "mercurius"
        sim.dt = 0.1313
        sim.ri_mercurius.safe_mode = 1
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)
        sim.save("test.bin")
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        sim = rebound.Simulation.from_file("test.bin")
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    
    

if __name__ == "__main__":
    unittest.main()
