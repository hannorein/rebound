import rebound
import unittest

class TestSerialize(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.t = 1.246
        self.sim.add(m=1.)
        self.sim.add(a=1.)
    
    def tearDown(self):
        self.sim = None

    def test_serialize(self):
        try:
            import numpy as np
        except:
            # Make numpy tests optional
            print("WARNING: Not testing serialization because numpy is not available")
            return;

        a = np.zeros((self.sim.N,3),dtype="float64")
        self.sim.serialize_particle_data(xyz=a)

        self.assertEqual(a[1][0],1)
        self.assertEqual(a[1][1],0)
        self.assertEqual(a[1][2],0)

        self.sim.particles[0].xyz = [0,0,0]

        b = np.zeros(self.sim.N,dtype="uint32")
        c = np.zeros(self.sim.N)
        self.sim.serialize_particle_data(r=c,hash=b)

        with self.assertRaises(AttributeError):
            self.sim.serialize_particle_data(r=b)

        with self.assertRaises(AttributeError):
            self.sim.serialize_particle_data(hash=a)
        
        with self.assertRaises(AttributeError):
            self.sim.serialize_particle_data(xyz=c)

    
    def test_set_serialize(self):
        try:
            import numpy as np
        except:
            # Make numpy tests optional
            print("WARNING: Not testing serialization because numpy is not available")
            return;
        
        for i in range(self.sim.N):
            self.sim.particles[i].xyz = [0,0,0]
            self.assertEqual(self.sim.particles[i].x,0)
            self.assertEqual(self.sim.particles[i].y,0)
            self.assertEqual(self.sim.particles[i].z,0)
        
        a = np.zeros((self.sim.N,3),dtype="float64")
        a[0][1] = 2
        a[1][2] = 6
        self.sim.set_serialized_particle_data(xyz=a)
        self.assertEqual(self.sim.particles[0].y,2)
        self.assertEqual(self.sim.particles[1].z,6)

    
if __name__ == "__main__":
    unittest.main()
