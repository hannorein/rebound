import rebound
import unittest
from ctypes import c_int
import math
import rebound.data
import warnings
    
    
class TestIntegratorWHFast512(unittest.TestCase):
    def test_whfast612_available(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value == 0:
            warnings.warn("WHFast512 not available. Cannot test WHFast512.", RuntimeWarning)
    def test_whfast512_com(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,9):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    sim = rebound.Simulation()
                    rebound.data.add_solar_system(sim,n_bodies)
                    sim.dt = sim.particles[1].P/23.456789
                    sim.move_to_com()
                    # Move to non-com frame.
                    for p in sim.particles:
                        p.vx += 0.1
                        p.y  += 0.1
                    sim.integrator = "whfast512"
                    sim.ri_whfast512.coordinates = coordinates
                    com0 = sim.com()
                    sim.integrate(1000.*2.*3.1415, exact_finish_time=0)
                    com1 = sim.com()
                    self.assertLess(math.fabs(com0.x+sim.t*com0.vx-com1.x),3.0e-13)
                    self.assertLess(math.fabs(com0.y+sim.t*com0.vy-com1.y),2.0e-16)
                    self.assertLess(math.fabs(com0.z+sim.t*com0.vz-com1.z),2.0e-16)
                    self.assertLess(math.fabs(com0.vx-com1.vx),1.0e-16)
                    self.assertLess(math.fabs(com0.vy-com1.vy),1.0e-16)
                    self.assertLess(math.fabs(com0.vz-com1.vz),1.0e-16)
    def test_whfast512_energy(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,9):
                sim = rebound.Simulation()
                rebound.data.add_solar_system(sim,n_bodies)
                sim.dt = sim.particles[1].P/23.456789
                sim.move_to_com()
                for p in sim.particles:
                    p.y  += 0.1
                    p.vx += 0.1
                sim.integrator = "whfast512"
                sim.ri_whfast512.coordinates = "jacobi"
                e0 = sim.energy()
                sim.integrate(1000.*2.*3.1415, exact_finish_time=0)
                e1 = sim.energy()
                self.assertLess(math.fabs((e0-e1)/e1),4.0e-9)
    
if __name__ == "__main__":
    unittest.main()
