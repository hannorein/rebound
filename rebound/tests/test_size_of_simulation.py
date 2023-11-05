import rebound
import unittest
from ctypes import c_size_t, sizeof

class TestSizeOfSimmulation(unittest.TestCase):
    
    def test_size_of_simulation(self):
        cl = rebound.clibrebound
        cl.reb_simulation_struct_size.res_type = c_size_t
        simulation_size_c = cl.reb_simulation_struct_size()

        self.assertEqual(simulation_size_c, sizeof(rebound.Simulation))

if __name__ == "__main__":
    unittest.main()
