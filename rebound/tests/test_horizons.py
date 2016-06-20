import rebound
import unittest
import datetime
import socket
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
            self.assertAlmostEqual(self.sim.particles[0].x,-0.17568959237103887,delta=1e-10)
            self.assertAlmostEqual(self.sim.particles[0].m,3.0404326480226416e-06,delta=1e-15)
        except socket.error: 
            print("Socket error. Most likely due to HORIZON being slow. Ignoring.")
            pass
    
    def test_notfound(self):
        with self.assertRaises(Exception):
            try:
                self.sim.add("BogusPlanet",date="2000-01-01 00:00")
            except socket.error: 
                print("Socket error. Most likely due to HORIZON being slow. Ignoring.")
                raise Exception("Socket error. Should have been bogus planet error. Ignoring")


if __name__ == "__main__":
    unittest.main()
