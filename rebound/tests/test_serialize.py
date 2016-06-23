import rebound
import unittest
import os
import math
import numpy as np

class TestSerialize(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.t = 1.246
        self.sim.add(m=1.)
        self.sim.add(a=1.)
    
    def tearDown(self):
        self.sim = None

    def test_serialize(self):
        a = np.zeros((self.sim.N,3),dtype="float64")
        self.sim.serialize_particle_data(xyz=a)

        self.assertEqual(a[1][0],1)
        self.assertEqual(a[1][1],0)
        self.assertEqual(a[1][2],0)

        b = np.zeros(self.sim.N,dtype="uint32")
        c = np.zeros(self.sim.N)
        self.sim.serialize_particle_data(r=c,hash=b)

        with self.assertRaises(AttributeError):
            self.sim.serialize_particle_data(r=b)

        with self.assertRaises(AttributeError):
            self.sim.serialize_particle_data(hash=a)
        
        with self.assertRaises(AttributeError):
            self.sim.serialize_particle_data(xyz=c)

    
if __name__ == "__main__":
    unittest.main()
