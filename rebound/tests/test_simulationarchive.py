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
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sim = rebound.Simulation.from_archive("test.bin")
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    def test_safe_mode_warning(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)
        x1 = sim.particles[1].x

        sa = rebound.SimulationArchive("test.bin")
        sim = sa.getSimulation(0.)
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.integrate(40.,exact_finish_time=0)
            x0 = sim.particles[1].x
            self.assertEqual(0, len(w)) # did not raise recalculating Jacobi warning

        self.assertEqual(x0,x1)
    
    def test_sa_whfasthelio_restart_safe_mode(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    def test_sa_whfasthelio_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "democraticheliocentric"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    def test_sa_whds_restart_safe_mode(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "whds"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "whds"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    def test_sa_whds_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "whds"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "whds"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    

    def test_sa_restart_safe_mode(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    def test_sa_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        t1 = sim.t
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x
        t0 = sim.t

        self.assertEqual(t0,t1)
        self.assertEqual(x0,x1)

        tget = 27.123
        sim = sa.getSimulation(tget,mode="exact");
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 25.123
        sim = sa.getSimulation(tget,mode="close");
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)

    def test_sa_restart_append(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.automateSimulationArchive("test.bin", 10.)
            sim.integrate(80.,exact_finish_time=0)
            self.assertEqual(1, len(w)) 

    
    def test_sa_restart_generator(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")

        times = [0.,11.,22.]
        for sim in sa.getSimulations(times,mode="close"):
            pass
        self.assertAlmostEqual(sim.t,22.058400000000105,delta=sim.dt)

    
    
    def test_sa_restart_corrector(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.ri_whfast.corrector = 5
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.ri_whfast.corrector = 5
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
        
        tget = 27.123
        sim = sa.getSimulation(tget,mode="exact");
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 25.123
        sim = sa.getSimulation(tget,mode="close");
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)
    
    
    def test_sa_restart_ias15(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
        
        tget = 27.123
        sim = sa.getSimulation(tget,mode="exact");
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 25.123
        sim = sa.getSimulation(tget,mode="close");
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)

    def test_sa_restart_ias15_walltime(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.automateSimulationArchive("test.bin", walltime = 0.01,deletefile=True)
        sim.integrate(1000.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        self.assertGreater(sim.t,100.)
        sim.integrate(2000.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.integrate(2000.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
        
        tget = 116.123
        sim = sa.getSimulation(tget,mode="exact");
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 116.123
        sim = sa.getSimulation(tget,mode="close");
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)


class TestSimulationArchiveWarningsErrors(unittest.TestCase):
    def test_sa_binary_missing(self):
        with self.assertRaises(ValueError):
            sa = rebound.SimulationArchive("testmissing.bin")
    def test_sa_binary_version(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.automateSimulationArchive("test.bin", walltime = 0.01,deletefile=True)
        sim.integrate(400.,exact_finish_time=0)
        with open("test.bin","r+b") as f:
            f.seek(30)
            f.write("1.0.0     ".encode('ascii'))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sa = rebound.SimulationArchive("test.bin")
            self.assertEqual(1, len(w)) 


class TestSimulationArchiveTmin(unittest.TestCase):
    def test_sa_tmin(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.integrate(400.,exact_finish_time=0)
        sim.automateSimulationArchive("test.bin", interval = 100.,deletefile=True)
        tmin = sim.t
        sim.integrate(800.,exact_finish_time=0)
        sa = rebound.SimulationArchive("test.bin")
        self.assertEqual(tmin,sa[0].t)
        self.assertEqual(tmin,sa.tmin)
        self.assertNotEqual(tmin,sa.tmax)
   

class TestSimulationArchiveMercurius(unittest.TestCase):
    def test_sa_mercurius_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "mercurius"
        sim.dt = 0.1313
        sim.ri_mercurius.safe_mode = 0
        sim.automateSimulationArchive("test.bin", 10.,deletefile=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "mercurius"
        sim.dt = 0.1313
        sim.ri_mercurius.safe_mode = 0
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

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

        sim = None
        sa = rebound.SimulationArchive("test.bin")
        sim = sa[-1]
        sim.integrate(80.,exact_finish_time=0)
        x1 = sim.particles[1].x
        
        
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "mercurius"
        sim.dt = 0.1313
        sim.ri_mercurius.safe_mode = 1
        sim.integrate(80.,exact_finish_time=0)
        x0 = sim.particles[1].x

        self.assertEqual(x0,x1)
    
    
    

if __name__ == "__main__":
    unittest.main()
