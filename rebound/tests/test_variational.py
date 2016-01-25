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
    def test_all_1st_order_init(self):
        vlist = ["a","e","i","Omega","omega","f","m"]
        for v in vlist:
            sim = rebound.Simulation()
            sim.add(m=1.)
            m=1e-3
            a=1.
            e=0.1
            inc=0.02
            Omega=0.3
            omega=0.56
            f=0.4
            sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
            sim.add(a=1.76)
            p = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
            vp = rebound.Particle(simulation=sim, primary=sim.particles[0],variation=v,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
            Delta = 1e-8
            if v=="a":
                a+=Delta
            if v=="e":
                e+=Delta
            if v=="i":
                inc+=Delta
            if v=="Omega":
                Omega+=Delta
            if v=="omega":
                omega+=Delta
            if v=="f":
                f+=Delta
            if v=="m":
                m+=Delta
            sp = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
            self.assertLess(abs((sp.x-p.x)/Delta-vp.x),1e-7)
            self.assertLess(abs((sp.y-p.y)/Delta-vp.y),1e-7)
            self.assertLess(abs((sp.z-p.z)/Delta-vp.z),1e-7)
            self.assertLess(abs((sp.vx-p.vx)/Delta-vp.vx),1e-7)
            self.assertLess(abs((sp.vy-p.vy)/Delta-vp.vy),1e-7)
            self.assertLess(abs((sp.vz-p.vz)/Delta-vp.vz),1e-7)
            self.assertLess(abs((sp.m-p.m)/Delta-vp.m),1e-7)
            
    def test_all_2nd_order_init(self):
        vlist = ["a","e","i","Omega","omega","f","m"]
        for v1 in vlist:
            for v2 in vlist:
                sim = rebound.Simulation()
                sim.add(m=1.)
                m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                sim.add(a=1.76)
                p = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                vp = rebound.Particle(simulation=sim, primary=sim.particles[0],variation=v1,variation2=v2,variation_order=2,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                Delta = 1e-4
                
                if v1==v2:

                    m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                    if v1=="a":
                        a+=Delta
                    if v1=="e":
                        e+=Delta
                    if v1=="i":
                        inc+=Delta
                    if v1=="Omega":
                        Omega+=Delta
                    if v1=="omega":
                        omega+=Delta
                    if v1=="f":
                        f+=Delta
                    if v1=="m":
                        m+=Delta
                    sp = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                    if v1=="a":
                        a-=Delta
                    if v1=="e":
                        e-=Delta
                    if v1=="i":
                        inc-=Delta
                    if v1=="Omega":
                        Omega-=Delta
                    if v1=="omega":
                        omega-=Delta
                    if v1=="f":
                        f-=Delta
                    if v1=="m":
                        m-=Delta
                    sm = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    prec = 1e-7
                    if v1=="m":
                        prec = 1e-2 # hard to find linear regime
                    self.assertLess(abs((sp.x-2.*p.x+sm.x)/Delta/Delta-vp.x),prec)
                    self.assertLess(abs((sp.y-2.*p.y+sm.y)/Delta/Delta-vp.y),prec)
                    self.assertLess(abs((sp.z-2.*p.z+sm.z)/Delta/Delta-vp.z),prec)
                    self.assertLess(abs((sp.vx-2.*p.vx+sm.vx)/Delta/Delta-vp.vx),prec)
                    self.assertLess(abs((sp.vy-2.*p.vy+sm.vy)/Delta/Delta-vp.vy),prec)
                    self.assertLess(abs((sp.vz-2.*p.vz+sm.vz)/Delta/Delta-vp.vz),prec)
                    self.assertLess(abs((sp.m-2.*p.m+sm.m)/Delta/Delta-vp.m),prec)
                else:
                    m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                    if v1=="a":
                        a+=Delta
                    if v1=="e":
                        e+=Delta
                    if v1=="i":
                        inc+=Delta
                    if v1=="Omega":
                        Omega+=Delta
                    if v1=="omega":
                        omega+=Delta
                    if v1=="f":
                        f+=Delta
                    if v1=="m":
                        m+=Delta
                    if v2=="a":
                        a+=Delta
                    if v2=="e":
                        e+=Delta
                    if v2=="i":
                        inc+=Delta
                    if v2=="Omega":
                        Omega+=Delta
                    if v2=="omega":
                        omega+=Delta
                    if v2=="f":
                        f+=Delta
                    if v2=="m":
                        m+=Delta
                    spp = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    
                    m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                    if v1=="a":
                        a+=Delta
                    if v1=="e":
                        e+=Delta
                    if v1=="i":
                        inc+=Delta
                    if v1=="Omega":
                        Omega+=Delta
                    if v1=="omega":
                        omega+=Delta
                    if v1=="f":
                        f+=Delta
                    if v1=="m":
                        m+=Delta
                    if v2=="a":
                        a-=Delta
                    if v2=="e":
                        e-=Delta
                    if v2=="i":
                        inc-=Delta
                    if v2=="Omega":
                        Omega-=Delta
                    if v2=="omega":
                        omega-=Delta
                    if v2=="f":
                        f-=Delta
                    if v2=="m":
                        m-=Delta
                    spm = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)


                    m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                    if v1=="a":
                        a-=Delta
                    if v1=="e":
                        e-=Delta
                    if v1=="i":
                        inc-=Delta
                    if v1=="Omega":
                        Omega-=Delta
                    if v1=="omega":
                        omega-=Delta
                    if v1=="f":
                        f-=Delta
                    if v1=="m":
                        m-=Delta
                    if v2=="a":
                        a+=Delta
                    if v2=="e":
                        e+=Delta
                    if v2=="i":
                        inc+=Delta
                    if v2=="Omega":
                        Omega+=Delta
                    if v2=="omega":
                        omega+=Delta
                    if v2=="f":
                        f+=Delta
                    if v2=="m":
                        m+=Delta
                    smp = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)

         
                    m,a,e,inc,Omega,omega,f=1e-3,1.,0.1,0.02,0.3,0.56,0.4
                    if v1=="a":
                        a-=Delta
                    if v1=="e":
                        e-=Delta
                    if v1=="i":
                        inc-=Delta
                    if v1=="Omega":
                        Omega-=Delta
                    if v1=="omega":
                        omega-=Delta
                    if v1=="f":
                        f-=Delta
                    if v1=="m":
                        m-=Delta
                    if v2=="a":
                        a-=Delta
                    if v2=="e":
                        e-=Delta
                    if v2=="i":
                        inc-=Delta
                    if v2=="Omega":
                        Omega-=Delta
                    if v2=="omega":
                        omega-=Delta
                    if v2=="f":
                        f-=Delta
                    if v2=="m":
                        m-=Delta
                    smm = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)

                    prec = 1e-7
                    if v1=="f" or v2=="f":
                        continue
                    if v1=="m" or v2=="m":
                        continue
                    print v1, v2
                    self.assertLess(abs((spp.x  -spm.x  -smp.x  +smm.x  )/Delta/Delta/4.-vp.x),prec)
                    self.assertLess(abs((spp.y  -spm.y  -smp.y  +smm.y  )/Delta/Delta/4.-vp.y),prec)
                    self.assertLess(abs((spp.z  -spm.z  -smp.z  +smm.z  )/Delta/Delta/4.-vp.z),prec)
                    self.assertLess(abs((spp.vx -spm.vx -smp.vx +smm.vx )/Delta/Delta/4.-vp.vx),prec)
                    self.assertLess(abs((spp.vy -spm.vy -smp.vy +smm.vy )/Delta/Delta/4.-vp.vy),prec)
                    self.assertLess(abs((spp.vz -spm.vz -smp.vz +smm.vz )/Delta/Delta/4.-vp.vz),prec)
                    self.assertLess(abs((spp.m  -spm.m  -smp.m  +smm.m  )/Delta/Delta/4.-vp.m),prec)
