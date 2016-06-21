import rebound
import unittest
import math
import numpy as np

class TestCollisions(unittest.TestCase):
    
    def test_tree_remove_both(self):
        sim = rebound.Simulation()
        boxsize = 50000.           
        sim.configure_box(boxsize)
        sim.integrator = "leapfrog"
        sim.boundary   = "open"
        sim.gravity    = "tree"
        sim.collision  = "tree"
        def cor_remove_both(r, c):
            r.contents.collisions_Nlog += 1
            return 3
        sim.collision_resolve = cor_remove_both
        
        while sim.N< 10:
            sim.add(m=1., r=100., x=np.random.uniform(-20,20),
                                y=np.random.uniform(-20,20),
                                z=np.random.uniform(-20,20))
        sim.dt = 0.001
        with self.assertRaises(rebound.NoParticles):
            sim.integrate(1000.)
        self.assertEqual(sim.collisions_Nlog,5)
    
    def test_direct_remove_both(self):
        sim = rebound.Simulation()
        boxsize = 50000.           
        sim.configure_box(boxsize)
        sim.integrator = "leapfrog"
        sim.boundary   = "open"
        sim.collision  = "direct"
        def cor_remove_both(r, c):
            r.contents.collisions_Nlog += 1
            return 3
        sim.collision_resolve = cor_remove_both
        
        while sim.N< 10:
            sim.add(m=1., r=100., x=np.random.uniform(-20,20),
                                y=np.random.uniform(-20,20),
                                z=np.random.uniform(-20,20))
        sim.dt = 0.001
        with self.assertRaises(rebound.NoParticles):
            sim.integrate(1000.)
        self.assertEqual(sim.collisions_Nlog,5)

    def test_tree_remove_one(self):
        sim = rebound.Simulation()
        boxsize = 50000.           
        sim.configure_box(boxsize)
        sim.integrator = "leapfrog"
        sim.boundary   = "open"
        sim.gravity    = "tree"
        sim.collision  = "tree"
        def cor_remove_both(r, c):
            r.contents.collisions_Nlog += 1
            return np.random.randint(1,3)
        sim.collision_resolve = cor_remove_both
        
        while sim.N< 50:
            sim.add(m=1., r=100., x=np.random.uniform(-2,2),
                                y=np.random.uniform(-2,2),
                                z=np.random.uniform(-2,2))
        sim.dt = 0.001
        sim.integrate(sim.dt)
        sim.integrate(2.*sim.dt)
        self.assertLess(sim.N,25)
    
    def test_direct_remove_one(self):
        sim = rebound.Simulation()
        boxsize = 50000.           
        sim.configure_box(boxsize)
        sim.integrator = "leapfrog"
        sim.boundary   = "open"
        sim.collision  = "direct"
        def cor_remove_both(r, c):
            r.contents.collisions_Nlog += 1
            return np.random.randint(1,3)
        sim.collision_resolve = cor_remove_both
        
        while sim.N< 50:
            sim.add(m=1., r=100., x=np.random.uniform(-2,2),
                                y=np.random.uniform(-2,2),
                                z=np.random.uniform(-2,2))
        sim.dt = 0.001
        sim.integrate(sim.dt)
        sim.integrate(2.*sim.dt)
        self.assertLess(sim.N,25)


if __name__ == "__main__":
    unittest.main()
