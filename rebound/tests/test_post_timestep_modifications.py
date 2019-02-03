import rebound
import unittest
import rebound.data as data

mdot = 1e-6
def ptm(sim):
    if sim.contents.ri_mercurius.mode == 0:
        sim.contents.particles[2].m -= sim.contents.dt_last_done*mdot

class TestPostTimestepModifications(unittest.TestCase):

    def test_ptm_ias15(self):
        sim = rebound.Simulation()
        sim.integrator = "ias15"
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        sim.post_timestep_modifications = ptm
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].m,1e-3-mdot*sim.t,delta=1e-13)

    def test_ptm_mercurius(self):
        sim = rebound.Simulation()
        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        sim.post_timestep_modifications = ptm
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].m,1e-3-mdot*sim.t,delta=1e-13)
    
    def test_ptm_mercurius_closeencounter(self):
        sim = rebound.Simulation()
        sim.integrator = "mercurius"
        sim.dt = 0.01
        sim.ri_mercurius.hillfac = 100000 # make sure encounter happens
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        sim.post_timestep_modifications = ptm
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].m,1e-3-mdot*sim.t,delta=1e-13)

    def test_ptm_whfast(self):
        sim = rebound.Simulation()
        sim.integrator = "whfast"
        sim.dt = 0.01
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        sim.post_timestep_modifications = ptm
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].m,1e-3-mdot*sim.t,delta=1e-13)

    def test_ptm_whfastjac(self):
        sim = rebound.Simulation()
        sim.integrator = "whfast"
        sim.ri_whfast.coordinates = "jacobi"
        sim.dt = 0.01
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=5)
        sim.move_to_com()
        sim.post_timestep_modifications = ptm
        sim.integrate(10.)
        self.assertAlmostEqual(sim.particles[2].m,1e-3-mdot*sim.t,delta=1e-13)


if __name__ == "__main__":
    unittest.main()

