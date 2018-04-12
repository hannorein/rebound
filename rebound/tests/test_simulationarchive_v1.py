import rebound
import unittest
import warnings

class TestSimulationArchive(unittest.TestCase):
    def test_sa_from_archive(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.simulationarchive_version = 1
        sim.initSimulationArchive("test.bin", 10.)
        sim.integrate(40.,exact_finish_time=0)
        print(sim.t)

        t1 = sim.t
        sim = None
        sim = rebound.Simulation.from_archive("test.bin")
        print(sim.t)
        sim.integrate(t1+1.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.integrate(t1+1.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    

if __name__ == "__main__":
    unittest.main()
