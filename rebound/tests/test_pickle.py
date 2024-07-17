import rebound
import unittest
import pickle
import warnings

class TestPickle(unittest.TestCase):
    def test_pickle_basic(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1,a=1)
        sim.step()
        
        with open("test.pickle", "wb") as f:
            pickle.dump(sim, f)

        with open("test.pickle", "rb") as f:
            sim2 = pickle.load(f)

        self.assertEqual(sim,sim2)
        sim.step()
        sim2.step()
        self.assertEqual(sim,sim2)
    
    def test_pickle_warning(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.collision_resolve = "merge"
        sim.step()
        
        with open("test.pickle", "wb") as f:
            pickle.dump(sim, f)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with open("test.pickle", "rb") as f:
                sim2 = pickle.load(f)

            self.assertEqual(1, len(w)) 
        self.assertNotEqual(sim,sim2)
        sim2.collision_resolve = "merge"
        sim.step()
        sim2.step()
        self.assertEqual(sim,sim2)

    def test_pickle_particle(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1,a=1)
        sim.step()
        
        with open("test.pickle", "wb") as f:
            pickle.dump(sim.particles[1], f)

        with open("test.pickle", "rb") as f:
            p2 = pickle.load(f)

        self.assertEqual(sim.particles[1], p2)

        # Pointers are set to zero when unpickling
        self.assertNotEqual(sim.particles[1]._sim, p2._sim)
        self.assertEqual(p2.sim, 0)

if __name__ == "__main__":
    unittest.main()
