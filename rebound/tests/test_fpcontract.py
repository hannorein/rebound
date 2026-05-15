import rebound
import unittest
from ctypes import c_int

class TestFPContract(unittest.TestCase):
    
    def test_fp_off(self):
        cl = rebound.clibrebound
        cl.reb_check_fp_contract.res_type = c_int
        fp_contract = cl.reb_check_fp_contract()

        # make sure floating point contraction are off
        self.assertEqual(fp_contract, 0)

class TestOpenMP(unittest.TestCase):
    
    def test_open_mp(self):
        # Will fail if not compiled with OpenMP enabled. Prints an error. Does not exit.
        rebound.omp_set_num_threads(10)

if __name__ == "__main__":
    unittest.main()
