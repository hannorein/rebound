import rebound
import unittest
import os
import warnings

class TestSimulationarchive(unittest.TestCase):
    def test_sa_step(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.save_to_file("test.bin", step=10,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sim = rebound.Simulation("test.bin")
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
    def test_sa_fromarchive(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 1
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sim = rebound.Simulation("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)
        x1 = sim.particles[1].x

        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim = sa.getSimulation(tget,mode="exact")
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 25.123
        sim = sa.getSimulation(tget,mode="close")
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)

    def test_sa_restart_append(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.ri_whfast.safe_mode = 0
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
        sim = sa[-1]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sim.save_to_file("test.bin", 10.)
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")

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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim = sa.getSimulation(tget,mode="exact")
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 25.123
        sim = sa.getSimulation(tget,mode="close")
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)
    
    
    def test_sa_restart_ias15(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim = sa.getSimulation(tget,mode="exact")
        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
        tget = 25.123
        sim = sa.getSimulation(tget,mode="close")
        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)

#    def test_sa_restart_ias15_walltime(self):
#        sim = rebound.Simulation()
#        sim.add(m=1)
#        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
#        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
#        sim.integrator = "ias15"
#        sim.dt = 0.1313
#        sim.save_to_file("test.bin", walltime = 0.1,delete_file=True)
#        sim.integrate(3000.,exact_finish_time=0)
#
#        sim = None
#        sa = rebound.Simulationarchive("test.bin")
#        sim = sa[-1]
#        self.assertGreater(sim.t,100.)
#        sim.integrate(20000.,exact_finish_time=0)
#        x1 = sim.particles[1].x
#        
#        
#        sim = rebound.Simulation()
#        sim.add(m=1)
#        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
#        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
#        sim.integrator = "ias15"
#        sim.dt = 0.1313
#        sim.integrate(20000.,exact_finish_time=0)
#        x0 = sim.particles[1].x
#
#        self.assertEqual(x0,x1)
#        
#        tget = 116.123
#        sim = sa.getSimulation(tget,mode="exact");
#        self.assertAlmostEqual(sim.t,tget,delta=1e-14)
#        tget = 116.123
#        sim = sa.getSimulation(tget,mode="close");
#        self.assertAlmostEqual(sim.t,tget,delta=sim.dt)
#

class TestSimulationarchiveWarningsErrors(unittest.TestCase):
    def test_sa_binary_missing(self):
        with self.assertRaises(RuntimeError):
            sa = rebound.Simulationarchive("testmissing.bin")
    def test_sa_binary_version(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "ias15"
        sim.dt = 0.1313
        sim.save_to_file("test.bin", walltime = 0.01,delete_file=True)
        sim.integrate(400.,exact_finish_time=0)
        with open("test.bin","r+b") as f:
            f.seek(30)
            f.write("1.0.0     ".encode('ascii'))
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sa = rebound.Simulationarchive("test.bin")
            self.assertEqual(1, len(w)) 


class TestSimulationarchiveTmin(unittest.TestCase):
    def test_sa_tmin(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "whfast"
        sim.dt = 0.1313
        sim.integrate(400.,exact_finish_time=0)
        sim.save_to_file("test.bin", interval = 100.,delete_file=True)
        tmin = sim.t
        sim.integrate(800.,exact_finish_time=0)
        sa = rebound.Simulationarchive("test.bin")
        self.assertEqual(tmin,sa[0].t)
        self.assertEqual(tmin,sa.tmin)
        self.assertNotEqual(tmin,sa.tmax)
   

class TestSimulationarchiveMercurius(unittest.TestCase):
    def test_sa_mercurius_restart(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        sim.integrator = "mercurius"
        sim.dt = 0.1313
        sim.ri_mercurius.safe_mode = 0
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(42.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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
        sim.save_to_file("test.bin", 10.,delete_file=True)
        sim.integrate(40.,exact_finish_time=0)

        sim = None
        sa = rebound.Simulationarchive("test.bin")
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

    def test_append_to_corrupt_snapshot(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3,a=1.)
        sim.add(m=5e-3,a=2.25)
    
        sim.save_to_file("simulationarchive.bin", interval=1000,delete_file=True) 
        sim.integrate(3001)
        with open('simulationarchive.bin', 'r+b') as f:
            f.seek(0, os.SEEK_END)          
            f.seek(f.tell() - 72, os.SEEK_SET)
            f.write(bytes(72)) # binary should be 15972 bytes, overwrite last 72 bytes with all zeros
        sa = rebound.Simulationarchive("simulationarchive.bin")
        sim = sa[-1]
        sim.save_to_file("simulationarchive.bin", interval=1000)
        sim.integrate(7001)
        sa = rebound.Simulationarchive("simulationarchive.bin")
        self.assertEqual(sa.nblobs, 8)
        self.assertAlmostEqual(sa[-1].t, 7000, places=0)
    

    def test_append_to_corrupt_snapshot_by_2bytes(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3,a=1.)
        sim.add(m=5e-3,a=2.25)
    
        sim.save_to_file("simulationarchive.bin", interval=1000,delete_file=True) 
        sim.integrate(3001)
        s1 = os.path.getsize('simulationarchive.bin')
        with open('simulationarchive.bin', 'r+b') as f:
            f.seek(0, os.SEEK_END)          
            f.seek(f.tell() - 2, os.SEEK_SET)
            f.truncate() # truncate by 2 bytes
        s2 = os.path.getsize('simulationarchive.bin')
        self.assertEqual(s1, s2+2)

        sa = rebound.Simulationarchive("simulationarchive.bin")
        sim = sa[-1]
        sim.save_to_file("simulationarchive.bin", interval=1000)
        sim.integrate(7001)
        sa = rebound.Simulationarchive("simulationarchive.bin")
        self.assertEqual(sa.nblobs, 8)
        self.assertAlmostEqual(sa[-1].t, 7000, places=0)
    
    def test_tree(self):
        sim = rebound.Simulation()
        sim.gravity = "tree"
        sim.integrator = "leapfrog"
        sim.configure_box(100,2,2,1)
        sim.add(m=1.)
        sim.add(m=1e-3,a=1.)
        sim.add(m=5e-3,a=2.25)
    
        sim.save_to_file("out.bin", interval=100,delete_file=True) 
        sim.integrate(305)
        
        sa = rebound.Simulationarchive("out.bin")
        self.assertEqual(sa.nblobs, 4)

        sim2 = sa[-1]
        sim2.integrate(305)
        self.assertEqual(sim, sim2)
    

if __name__ == "__main__":
    unittest.main()
