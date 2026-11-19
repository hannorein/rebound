import rebound
import unittest
import os
import warnings

class TestHorizons(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
    
    def tearDown(self):
        self.sim = None
    
    def test_earth(self):
        try:
            with warnings.catch_warnings(record=True) as w: 
                warnings.simplefilter("always")
                self.sim.add("Earth",date="2000-01-01 00:00")
                self.assertEqual(0,len(w))
            self.assertAlmostEqual(self.sim.particles[0].x,-0.17569031580176828,delta=1e-10)
            self.assertAlmostEqual(self.sim.particles[0].m,3.0404326480226416e-06,delta=1e-15)
        except: 
            # Output for GitHub actions:
            if os.getenv('GITHUB_ACTIONS') == 'true':
                print("::warning file=rebound/tests/test_horizons.py,line=24:: HORIZONS error. Most likely due to HORIZON being slow.")
            else:
                raise Exception("HORIZONS error. Most likely due to HORIZON being slow.")

if __name__ == "__main__":
    unittest.main()
