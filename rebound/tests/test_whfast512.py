import rebound
import unittest
from ctypes import c_int
import math
import rebound.data
import warnings
    
def gr_force(simp):
    C2 = 10065.32 * 10065.32;
    # Inefficient implementation of the GR force for WHFast.
    source = simp.contents.particles[0]
    prefac1 = 6.*(simp.contents.G*source.m)*(simp.contents.G*source.m)/C2
    for i in range(1,simp.contents.N):
        dx = simp.contents.particles[i].x - source.x
        dy = simp.contents.particles[i].y - source.y
        dz = simp.contents.particles[i].z - source.z
        r2 = dx*dx + dy*dy + dz*dz
        prefac = prefac1/(r2*r2)

        simp.contents.particles[i].ax -= prefac*dx
        simp.contents.particles[i].ay -= prefac*dy
        simp.contents.particles[i].az -= prefac*dz
        simp.contents.particles[0].ax += simp.contents.particles[i].m/source.m*prefac*dx
        simp.contents.particles[0].ay += simp.contents.particles[i].m/source.m*prefac*dy
        simp.contents.particles[0].az += simp.contents.particles[i].m/source.m*prefac*dz
    
class TestIntegratorWHFast512(unittest.TestCase):
    def test_whfast612_available(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value == 0:
            warnings.warn("WHFast512 not available. Cannot test WHFast512. Try `export AVX512=1`, then `pip install`.", RuntimeWarning)

    def test_whfast512_gr(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,10):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    sim = rebound.Simulation()
                    rebound.data.add_solar_system(sim,n_bodies)
                    sim.dt = sim.particles[1].P/23.456789
                    sim.move_to_com()
                    # Move to non-com frame to make the test a bit harder.
                    for p in sim.particles:
                        p.vx += 0.1
                        p.y  += 0.1
                    sim2 = sim.copy()

                    sim.integrator = "whfast512"
                    sim.ri_whfast512.coordinates = coordinates
                    sim.ri_whfast512.gr_potential = 1
                    sim.exact_finish_time = 0
                    sim.steps(10)
                    sim.synchronize()
                    
                    sim2.integrator = "whfast"
                    sim2.ri_whfast.coordinates = coordinates
                    sim2.ri_whfast.safe_mode = 0
                    sim2.additional_forces = gr_force
                    sim2.steps(10)
                    sim2.synchronize()
                    
                    for p1, p2 in zip(sim.particles, sim2.particles):
                        self.assertLess(math.fabs(p1.x-p2.x),8e-13)
                        self.assertLess(math.fabs(p1.y-p2.y),8e-13)
                        self.assertLess(math.fabs(p1.z-p2.z),8e-13)
                        self.assertLess(math.fabs(p1.vx-p2.vx),8e-13)
                        self.assertLess(math.fabs(p1.vy-p2.vy),8e-13)
                        self.assertLess(math.fabs(p1.vz-p2.vz),8e-13)

    def test_whfast512_com(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,10):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    sim = rebound.Simulation()
                    rebound.data.add_solar_system(sim,n_bodies)
                    sim.dt = sim.particles[1].P/23.456789
                    sim.move_to_com()
                    # Move to non-com frame to make the test a bit harder.
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
            for n_bodies in range(2,10):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    sim = rebound.Simulation()
                    rebound.data.add_solar_system(sim,n_bodies)
                    sim.dt = sim.particles[1].P/23.456789
                    sim.move_to_com()
                    # Move to non-com frame to make the test a bit harder.
                    for p in sim.particles:
                        p.y  += 0.1
                        p.vx += 0.1
                    sim.integrator = "whfast512"
                    sim.ri_whfast512.coordinates = coordinates
                    e0 = sim.energy()
                    sim.integrate(1000.*2.*3.1415, exact_finish_time=0)
                    e1 = sim.energy()
                    self.assertLess(math.fabs((e0-e1)/e1),4.0e-9)
    def test_whfast512_whfast(self):
        # WHfast512 should give same results as WHFast
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,10):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    sim = rebound.Simulation()
                    rebound.data.add_solar_system(sim,n_bodies)
                    sim.dt = sim.particles[1].P/23.456789
                    sim.move_to_com()
                    # Move to non-com frame to make the test a bit harder.
                    for p in sim.particles:
                        p.y  += 0.1
                        p.vx += 0.1
                    sim2 = sim.copy()

                    sim.integrator = "whfast512"
                    sim.ri_whfast512.coordinates = coordinates
                    sim.exact_finish_time = 0
                    sim.steps(10)
                    sim.synchronize()
                    
                    sim2.integrator = "whfast"
                    sim2.ri_whfast.coordinates = coordinates
                    sim2.ri_whfast.safe_mode = 0
                    sim2.steps(10)
                    sim2.synchronize()
                    for p1, p2 in zip(sim.particles, sim2.particles):
                        self.assertLess(math.fabs(p1.x-p2.x),8e-15)
                        self.assertLess(math.fabs(p1.y-p2.y),8e-15)
                        self.assertLess(math.fabs(p1.z-p2.z),8e-15)
                        self.assertLess(math.fabs(p1.vx-p2.vx),8e-15)
                        self.assertLess(math.fabs(p1.vy-p2.vy),8e-15)
                        self.assertLess(math.fabs(p1.vz-p2.vz),8e-15)
    def test_whfast512_restart(self):
        # Exact reproducibility
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,10):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    for keep_unsynchronized in [1, 0]:
                        sim = rebound.Simulation()
                        rebound.data.add_solar_system(sim,n_bodies)
                        sim.dt = sim.particles[1].P/23.456789
                        sim.move_to_com()
                        # Move to non-com frame to make the test a bit harder.
                        for p in sim.particles:
                            p.y  += 0.1
                            p.vx += 0.1

                        sim.integrator = "whfast512"
                        sim.ri_whfast512.coordinates = coordinates
                        sim.ri_whfast512.keep_unsynchronized = keep_unsynchronized
                        sim.exact_finish_time = 0
                        sim.steps(10)
                        sim.save_to_file("test1.bin",delete_file=True)
                        sim.synchronize()
                        sim.save_to_file("test2.bin",delete_file=True)
                        sim.steps(10)
                        sim.synchronize()

                        sim1 = rebound.Simulation("test1.bin")
                        sim1.synchronize()
                        sim1.steps(10)
                        sim1.synchronize()
                        self.assertEqual(sim, sim1)
                       
                        sim2 = rebound.Simulation("test2.bin")
                        sim2.steps(10)
                        sim2.synchronize()
                        self.assertEqual(sim, sim2)
    def test_whfast512_keep_unsynchronized(self):
        if c_int.in_dll(rebound.clibrebound, "reb_integrator_whfast512_available").value:
            for n_bodies in range(2,10):
                for coordinates in ["jacobi", "democraticheliocentric"]:
                    for keep_unsynchronized in [1, 0]:
                        sim = rebound.Simulation()
                        rebound.data.add_solar_system(sim,n_bodies)
                        sim.dt = sim.particles[1].P/23.456789
                        sim.move_to_com()
                        # Move to non-com frame to make the test a bit harder.
                        for p in sim.particles:
                            p.y  += 0.1
                            p.vx += 0.1

                        sim.integrator = "whfast512"
                        sim.ri_whfast512.coordinates = coordinates
                        sim.ri_whfast512.keep_unsynchronized = keep_unsynchronized
                        sim.exact_finish_time = 0
                        sim.steps(10)
                        sim.save_to_file("test1.bin",delete_file=True)
                        sim.synchronize() # extra synchronize not in sim1
                        sim.steps(10)
                        sim.synchronize()

                        sim1 = rebound.Simulation("test1.bin")
                        sim1.steps(10)
                        sim1.synchronize()
                        if keep_unsynchronized:
                            self.assertEqual(sim, sim1)
                        else:
                            self.assertNotEqual(sim, sim1)
                   
    
if __name__ == "__main__":
    unittest.main()
