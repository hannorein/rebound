import rebound
import unittest
import datetime

class TestVariationalIAS15(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3,a=1)
        self.sim.add(a=1.76)
        self.sim.move_to_com()
        
        self.DeltaX = 0.001
        
        simshifted = rebound.Simulation()
        simshifted.add(m=1.)
        simshifted.add(m=1.e-3,a=1)
        simshifted.add(a=1.76)
        simshifted.particles[1].x += self.DeltaX
        simshifted.move_to_com()
        simshifted.integrate(100.)
        self.Xshifted = simshifted.particles[2].x
    
    def tearDown(self):
        self.sim = None
    
    def test_var_1st_order(self):
        var_i = self.sim.add_variational()
        self.sim.particles[var_i+1].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*self.sim.particles[var_i+2].x 
        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-3)
    
    def test_var_2nd_order(self):
        var_i = self.sim.add_variational()
        var_ii = self.sim.add_variational(order=2,index_1st_order_a=var_i, index_1st_order_b=var_i)
        self.sim.particles[var_i+1].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*self.sim.particles[var_i+2].x + self.DeltaX**2/2.*self.sim.particles[var_ii+2].x 

        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-5)
    

class TestVariationalWHFast(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.integrator = "WHFast"
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3,a=1)
        self.sim.add(a=1.76)
        self.sim.move_to_com()
        
        self.DeltaX = 0.001
        
        simshifted = rebound.Simulation()
        simshifted.integrator = "WHFast"
        simshifted.add(m=1.)
        simshifted.add(m=1.e-3,a=1)
        simshifted.add(a=1.76)
        simshifted.particles[1].x += self.DeltaX
        simshifted.move_to_com()
        simshifted.integrate(100.)
        self.Xshifted = simshifted.particles[2].x
    
    def tearDown(self):
        self.sim = None
    
    def test_var_1st_order(self):
        var_i = self.sim.add_variational()
        self.sim.particles[var_i+1].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*self.sim.particles[var_i+2].x 
        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-3)
    
    


class TestVariationalTestparticleIAS15(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3,a=1)
        self.sim.add(a=1.76)
        self.sim.move_to_com()
        
        self.DeltaX = 0.001
        
        simshifted = rebound.Simulation()
        simshifted.add(m=1.)
        simshifted.add(m=1.e-3,a=1)
        simshifted.add(a=1.76)
        simshifted.particles[2].x += self.DeltaX
        simshifted.move_to_com()
        simshifted.integrate(100.)
        self.Xshifted = simshifted.particles[2].x
    
    def tearDown(self):
        self.sim = None
    
    def test_var_1st_order(self):
        var_i = self.sim.add_variational(testparticle=2)
        self.sim.particles[var_i].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*self.sim.particles[var_i].x 
        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-2)
    
    def test_var_2nd_order(self):
        var_i = self.sim.add_variational(testparticle=2)
        var_ii = self.sim.add_variational(order=2,index_1st_order_a=var_i, index_1st_order_b=var_i, testparticle=2)
        self.sim.particles[var_i].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*self.sim.particles[var_i].x + self.DeltaX**2/2.*self.sim.particles[var_ii].x 

        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-3)
    

