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

if __name__ == "__main__":
    unittest.main()
