import rebound
import unittest
import rebound.data as data

class TestAdditionalForces(unittest.TestCase):
    def test_af_ias15(self):
        sim = rebound.Simulation()
        sim.integrator = "ias15"
        sim.force_is_velocity_dependent = 1
        sim.add(m=1)
        sim.add(m=1e-6,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        def af(sim):
            fac = 0.01
            sim.contents.particles[2].ax -= fac*sim.contents.particles[2].vx 
            sim.contents.particles[2].ay -= fac*sim.contents.particles[2].vz 
            sim.contents.particles[2].az -= fac*sim.contents.particles[2].vy 
        sim.additional_forces = af
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].a,4.86583,delta=1e-5)

    def test_af_mercurius(self):
        sim = rebound.Simulation()
        sim.integrator = "mercurius"
        sim.dt = 0.005
        sim.force_is_velocity_dependent = 1
        sim.add(m=1)
        sim.add(m=1e-6,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        def af(sim):
            fac = 0.01
            sim.contents.particles[2].ax -= fac*sim.contents.particles[2].vx 
            sim.contents.particles[2].ay -= fac*sim.contents.particles[2].vz 
            sim.contents.particles[2].az -= fac*sim.contents.particles[2].vy 
        sim.additional_forces = af
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].a,4.86583,delta=1e-4)
    
    def test_af_mercurius_closeencounter(self):
        sim = rebound.Simulation()
        sim.integrator = "mercurius"
        sim.ri_mercurius.hillfac = 100000 # make sure encounter happens
        sim.dt = 0.005
        sim.force_is_velocity_dependent = 1
        sim.add(m=1)
        sim.add(m=1e-6,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        def af(sim):
            fac = 0.01
            sim.contents.particles[2].ax -= fac*sim.contents.particles[2].vx 
            sim.contents.particles[2].ay -= fac*sim.contents.particles[2].vz 
            sim.contents.particles[2].az -= fac*sim.contents.particles[2].vy 
        sim.additional_forces = af
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].a,4.86583,delta=1e-4)


    def test_af_whfast(self):
        sim = rebound.Simulation()
        sim.integrator = "whfast"
        sim.dt = 0.01
        sim.force_is_velocity_dependent = 1
        sim.add(m=1)
        sim.add(m=1e-6,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        def af(sim):
            fac = 0.01
            sim.contents.particles[2].ax -= fac*sim.contents.particles[2].vx 
            sim.contents.particles[2].ay -= fac*sim.contents.particles[2].vz 
            sim.contents.particles[2].az -= fac*sim.contents.particles[2].vy 
        sim.additional_forces = af
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].a,4.86583,delta=1e-4)

    def test_af_whfastjac(self):
        sim = rebound.Simulation()
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "jacobi"
        sim.dt = 0.01
        sim.force_is_velocity_dependent = 1
        sim.add(m=1)
        sim.add(m=1e-6,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        def af(sim):
            fac = 0.01
            sim.contents.particles[2].ax -= fac*sim.contents.particles[2].vx 
            sim.contents.particles[2].ay -= fac*sim.contents.particles[2].vz 
            sim.contents.particles[2].az -= fac*sim.contents.particles[2].vy 
        sim.additional_forces = af
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].a,4.86583,delta=1e-4)


if __name__ == "__main__":
    unittest.main()

