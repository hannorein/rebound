import rebound
import unittest
import math
import numpy as np

class TestLineCollisions(unittest.TestCase):
    
    def test_direct_miss(self):
        # Should miss the collision
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.collision  = "direct"
        sim.dt = 10
        sim.add(r=1,x=0)
        sim.add(r=1,x=3,vx=-1)
        sim.integrate(10)
    def test_line_find(self):
        # Should find the collision
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.collision  = "line"
        sim.dt = 10
        sim.add(r=1,x=0)
        sim.add(r=1,x=3,vx=-1)
        with self.assertRaises(rebound.Collision) as context:
            sim.integrate(10)
    def test_line_miss_overlap(self):
        # Should miss the collision because overlapping at t=0
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.collision  = "line"
        sim.dt = 10
        sim.add(r=1,x=0)
        sim.add(r=1,x=1,vx=-1)
        sim.integrate(10)
    def test_line_find_overlap(self):
        # Should find the collision because not overlapping at t=0, only at end
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.collision  = "line"
        sim.dt = 10
        sim.add(r=1,x=0)
        sim.add(r=1,x=11,vx=-1)
        with self.assertRaises(rebound.Collision) as context:
            sim.integrate(10)


class TestCollisions(unittest.TestCase):
    
    def test_tree_remove_both(self):
        sim = rebound.Simulation()
        boxsize = 50000.           
        sim.configure_box(boxsize)
        sim.integrator = "leapfrog"
        sim.boundary   = "open"
        sim.gravity    = "tree"
        sim.collision  = "tree"
        sim.collision_resolve = "hardsphere"
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
