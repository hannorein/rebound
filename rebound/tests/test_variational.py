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
    

class TestVariationalInitTests(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1.)
        self.sim.add(m=1.e-3,a=1,e=0.1,inc=0.02,Omega=0.3,omega=0.56,f=0.4)
        self.sim.add(a=1.76)
    def test_all_1st_order_init(self):
        vlist = ["a","e","i","Omega","omega","f","m"]
        for v in vlist:
            p = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1,e=0.1,inc=0.02,Omega=0.3,omega=0.56,f=0.4)
            vp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],variation=v,m=1.e-3,a=1,e=0.1,inc=0.02,Omega=0.3,omega=0.56,f=0.4)
            Delta = 1e-8
            if v=="a":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1+Delta,e=0.1,inc=0.02,Omega=0.3,omega=0.56,f=0.4)
            if v=="e":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1,e=0.1+Delta,inc=0.02,Omega=0.3,omega=0.56,f=0.4)
            if v=="i":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1,e=0.1,inc=0.02+Delta,Omega=0.3,omega=0.56,f=0.4)
            if v=="Omega":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1,e=0.1,inc=0.02,Omega=0.3+Delta,omega=0.56,f=0.4)
            if v=="omega":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1,e=0.1,inc=0.02,Omega=0.3,omega=0.56+Delta,f=0.4)
            if v=="f":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3,a=1,e=0.1,inc=0.02,Omega=0.3,omega=0.56,f=0.4+Delta)
            if v=="m":
                sp = rebound.Particle(simulation=self.sim, primary=self.sim.particles[0],m=1.e-3+Delta,a=1,e=0.1,inc=0.02,Omega=0.3,omega=0.56,f=0.4)
            self.assertLess(abs((sp.x-p.x)/Delta-vp.x),1e-7)
            self.assertLess(abs((sp.y-p.y)/Delta-vp.y),1e-7)
            self.assertLess(abs((sp.z-p.z)/Delta-vp.z),1e-7)
            self.assertLess(abs((sp.vx-p.vx)/Delta-vp.vx),1e-7)
            self.assertLess(abs((sp.vy-p.vy)/Delta-vp.vy),1e-7)
            self.assertLess(abs((sp.vz-p.vz)/Delta-vp.vz),1e-7)
            self.assertLess(abs((sp.m-p.m)/Delta-vp.m),1e-7)
            
    def tearDown(self):
        self.sim = None
