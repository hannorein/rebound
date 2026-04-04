import rebound
import unittest
from ctypes import c_size_t, sizeof, c_uint32, POINTER, c_char_p
from rebound.simulation import Integrator

class TestSizeOfSimulation(unittest.TestCase):
    
    def test_size_of_simulation(self):
        cl = rebound.clibrebound
        cl.reb_simulation_struct_size.res_type = c_size_t
        simulation_size_c = cl.reb_simulation_struct_size()

        self.assertEqual(simulation_size_c, sizeof(rebound.Simulation))

class TestUniqueIntegrator(unittest.TestCase):
    def setUp(self):
        clibrebound = rebound.clibrebound
        integrators_available_N = c_uint32.in_dll(clibrebound, "reb_integrators_available_N").value
        integrators_available = (POINTER(Integrator) * integrators_available_N).in_dll(clibrebound, "reb_integrators_available")
        integrators_available_names = (c_char_p * integrators_available_N).in_dll(clibrebound, "reb_integrators_available_names")
        self.INTEGRATORS = dict(zip([n.decode("utf-8") for n in integrators_available_names], [i.contents for i in integrators_available]))
    def test_unique_integrator_names(self):
        self.assertEqual(len(self.INTEGRATORS.keys()), len(set(self.INTEGRATORS.keys())))
    def test_unique_integrator_ids(self):
        self.assertEqual(len(self.INTEGRATORS.values()), len(set([i.id for i in self.INTEGRATORS.values()])))


if __name__ == "__main__":
    unittest.main()
