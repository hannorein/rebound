import rebound
import unittest
import math

class TestIntegrator(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.move_to_com()
    
    def tearDown(self):
        self.sim = None

    def test_force_exact_finish_time(self):
        self.sim.integrator = "modleapfrog"
        self.sim.integrator.order = 2
        self.sim.dt = 1
        self.sim.steps(1)
        exact_finish_time = self.sim.exact_finish_time
        self.assertEqual(exact_finish_time,0)
    
    def test_modleapfrog_wrong_order(self):
        self.sim.integrator = "modleapfrog"
        self.sim.integrator.order = 42
        with self.assertRaises(RuntimeError):
            self.sim.steps(1)

    def test_modleapfrog_energy_solar_system_order(self):
        for order, target in zip([2,4,6,8], [1e-5, 2e-6, 1e-9, 1e-12]):
            self.sim.integrator = "modleapfrog"
            self.sim.integrator.order = order
            self.sim.add("solar system")
            self.sim.dt = 1e-5
            e0 = self.sim.energy()
            self.assertNotEqual(e0,0.)
            self.sim.integrate(1e3)
            e1 = self.sim.energy()
            self.assertLess(math.fabs((e0-e1)/e1), target)
    
    def test_modleapfrog_energy_3body(self):
        for order, target in zip([2,4,6,8], [1e-4, 2e-7, 1e-11, 1e-11]):
            self.sim.integrator = "modleapfrog"
            self.sim.integrator.order = order
            self.sim.dt = 1e-2

            primary = rebound.Particle(m=1000)
            self.sim.add(primary)
            secondary = rebound.Particle(simulation=self.sim, primary=primary, m=0.1, a =1, e=0.999)
            tertiary = rebound.Particle(simulation=self.sim, primary=primary, m=10.0, a =10, e=0.2)
            self.sim.add(secondary)
            self.sim.add(tertiary)
            orbits = self.sim.orbits()
            T = 1.3* orbits[0].P

            e0 = self.sim.energy()
            self.assertNotEqual(e0,0.)
            self.sim.integrate(T)
            e1 = self.sim.energy()
            self.assertLess(math.fabs((e0-e1)/e1), target)

    def test_modleapfrog_energy_collision(self):
        self.sim.integrator = "modleapfrog"
        self.sim.integrator.order = 4
        self.sim.dt = 1e-1

        self.sim.add(rebound.Particle(m=1.0, x=1))
        self.sim.add(rebound.Particle(m=1.0, x=-1))

        e0 = self.sim.energy()
        self.assertNotEqual(e0,0.)
        self.sim.integrate(1e3)
        e1 = self.sim.energy()
        self.assertLess(math.fabs((e0-e1)/e1), 1e-13)


    def test_modleapfrog_dt_last_done(self):
        self.sim.integrator = "modleapfrog"
        self.sim.integrator.order = 4

        dtau = 1e-1
        self.sim.dt = dtau

        self.sim.add(rebound.Particle(m=1.0, x=1))
        self.sim.add(rebound.Particle(m=1.0, x=-1))

        self.sim.steps(1)

        self.assertNotEqual(self.sim.dt_last_done, dtau)
        
if __name__ == "__main__":
    unittest.main()
