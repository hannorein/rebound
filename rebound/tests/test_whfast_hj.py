import ctypes
import unittest

import rebound


class TestIntegratorWHFastHJGivenTree(unittest.TestCase):
    def make_sim(self):
        sim = rebound.Simulation()
        sim.G = 1.0
        sim.integrator = "whfast_hj"
        sim.dt = 0.01
        sim.add(m=1.0)
        sim.add(m=0.1, a=1.0, e=0.05)
        sim.add(m=0.0, a=2.0, e=0.1)
        sim.move_to_com()
        return sim

    def hj_tree_string(self, sim):
        clibrebound = rebound.clibrebound
        clibrebound.reb_integrator_whfast_hj_tree_to_string.argtypes = [
            ctypes.POINTER(rebound.Simulation),
            ctypes.c_char_p,
            ctypes.c_size_t,
        ]
        clibrebound.reb_integrator_whfast_hj_tree_to_string.restype = ctypes.c_int

        small = ctypes.create_string_buffer(1)
        required = clibrebound.reb_integrator_whfast_hj_tree_to_string(
            ctypes.byref(sim),
            small,
            ctypes.c_size_t(1),
        )
        self.assertGreaterEqual(required, 0)

        buffer = ctypes.create_string_buffer(required + 1)
        clibrebound.reb_integrator_whfast_hj_tree_to_string(
            ctypes.byref(sim),
            buffer,
            ctypes.c_size_t(required + 1),
        )
        return buffer.value.decode("ascii")

    def test_given_tree_matches_automatic_tree(self):
        sim_auto = self.make_sim()
        sim_given = self.make_sim()

        tree = self.hj_tree_string(sim_auto)
        sim_auto.integrate(0.1, exact_finish_time=0)
        sim_given.integrate(0.1, exact_finish_time=0, given_tree=True, tree=tree)

        self.assertEqual(sim_given.integrator.given_tree, 1)
        for particle_auto, particle_given in zip(sim_auto.particles, sim_given.particles):
            for attr in ("x", "y", "z", "vx", "vy", "vz"):
                self.assertAlmostEqual(
                    getattr(particle_auto, attr),
                    getattr(particle_given, attr),
                    delta=1e-14,
                )

    def test_given_tree_accepts_nested_pairs(self):
        sim = self.make_sim()
        sim.integrate(0.01, exact_finish_time=0, given_tree=True, tree=[[1, 2], 3])
        self.assertEqual(sim.integrator.given_tree, 1)

    def test_given_tree_accepts_binary_plus_particles_mode(self):
        sim_auto = self.make_sim()
        sim_mode = self.make_sim()

        tree = self.hj_tree_string(sim_auto)
        sim_auto.integrate(0.1, exact_finish_time=0, given_tree=True, tree=tree)
        sim_mode.integrate(0.1, exact_finish_time=0, given_tree=True, tree="binary_plus_particles")

        self.assertEqual(sim_mode.integrator.given_tree, 1)
        for particle_auto, particle_mode in zip(sim_auto.particles, sim_mode.particles):
            for attr in ("x", "y", "z", "vx", "vy", "vz"):
                self.assertAlmostEqual(
                    getattr(particle_auto, attr),
                    getattr(particle_mode, attr),
                    delta=1e-14,
                )

    def test_given_tree_rejects_duplicate_particle(self):
        sim = self.make_sim()
        with self.assertRaises(RuntimeError):
            sim.integrate(0.01, exact_finish_time=0, given_tree=True, tree="[1,1]")


if __name__ == "__main__":
    unittest.main()
