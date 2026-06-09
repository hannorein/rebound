import rebound
import unittest
import math
import warnings
    
def gr_potential(sim):
    C2 = 10065.32 * 10065.32
    source = sim.particles[0]
    mu = sim.G*source.m;
    prefac = 3.*mu*mu/C2
    H = 0.0

    for p in sim.particles[1:]:
        dx = p.x - source.x
        dy = p.y - source.y
        dz = p.z - source.z
        r2 = dx*dx + dy*dy + dz*dz
        H -= prefac*p.m/r2
    return H
    
class TestIntegratorWHFast512(unittest.TestCase):
    def test_whfast512_basic(self):
        if not rebound.avx512_available: return
        sim = rebound.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy()
        sim.steps(1)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e0),5e-11)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy()
        self.assertEqual(sim.t, sim.dt*(1+1e5))
        self.assertLess(math.fabs((e0-e1)/e0),2e-9)

    def test_whfast512_corrector(self):
        if not rebound.avx512_available: return
        sim = rebound.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.corrector = 17
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy()
        sim.steps(1)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e0),6e-15)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e0),2e-12)

    def test_whfast512_gr(self):
        if not rebound.avx512_available: return
        sim = rebound.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.corrector = 17
        sim.integrator.gr_potential = 1
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy() + gr_potential(sim)
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),6e-15)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),2e-12)

    
    def test_whfast512_fewer_than_8_particles(self):
        if not rebound.avx512_available: return
        sim = rebound.Simulation()
        sim.add("outer solar system")
        sim.integrator = "whfast512"
        sim.integrator.corrector = 17
        sim.integrator.gr_potential = 1
        sim.dt = 20.0/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy() + gr_potential(sim)
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),4e-14)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),1e-11)


if __name__ == "__main__":
    unittest.main()
