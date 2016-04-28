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
        var_i = self.sim.add_variation()
        var_i.particles[1].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*var_i.particles[2].x 
        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-3)
    
    def test_var_2nd_order(self):
        var_i = self.sim.add_variation()
        var_ii = self.sim.add_variation(order=2,first_order=var_i)
        var_i.particles[1].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*var_i.particles[2].x + self.DeltaX**2/2.*var_ii.particles[2].x 

        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-5)
    
    def test_var_restart(self):
        var_i = self.sim.add_variation()
        var_ii = self.sim.add_variation(order=2,first_order=var_i)
        var_i.particles[1].x = 1.
        self.sim.save("test.bin")
        self.sim = None
        sim = rebound.Simulation.from_file("test.bin")
        var_i = sim.var_config[0]
        var_ii = sim.var_config[1]
        var_i.particles[1].x = 1.
        sim.integrate(100.)
        X_var = sim.particles[2].x + self.DeltaX*var_i.particles[2].x + self.DeltaX**2/2.*var_ii.particles[2].x 

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
        var_i = self.sim.add_variation()
        var_i.particles[1].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*var_i.particles[2].x 
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
        var_i = self.sim.add_variation(testparticle=2)
        var_i.particles[0].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*var_i.particles[0].x 
        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-2)
    
    def test_var_2nd_order(self):
        var_i = self.sim.add_variation(testparticle=2)
        var_ii = self.sim.add_variation(order=2,first_order=var_i, testparticle=2)
        var_i.particles[0].x = 1.
        self.sim.integrate(100.)
        X_var = self.sim.particles[2].x + self.DeltaX*var_i.particles[0].x + self.DeltaX**2/2.*var_ii.particles[0].x 

        self.assertAlmostEqual(self.Xshifted,X_var,delta=1e-3)
    

class TestVariationalInitTests(unittest.TestCase):
    def test_all_1st_order_init(self):
        vlist = ["a","e","i","Omega","omega","f","m"]
        paramslist = [ 
                (1e-3,1.,0.1,0.02,0.3,0.56,0.4),
                (1e-6,2.,0.02,0.0132,0.33,1.56,0.14),
                (234.3e-6,1.7567,0.561,0.572,0.573,2.56,0.354),
                (1e-2,1.7567,0.1561,0.15472,0.24573,12.56,1.354),
                (1e-7,3.7567,0.00061,0.23572,0.523473,2.56,3.354),
                ]
        for params in paramslist:
            for v in vlist:
                sim = rebound.Simulation()
                sim.add(m=1.)
                m,a,e,inc,Omega,omega,f= params
                Delta=1e-8
                sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                sim.add(a=1.76)
                p = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                vp = rebound.Particle(simulation=sim, primary=sim.particles[0],variation=v,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
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
                self.assertLess(abs((sp.x-p.x)/Delta-vp.x),1e-6)
                self.assertLess(abs((sp.y-p.y)/Delta-vp.y),1e-6)
                self.assertLess(abs((sp.z-p.z)/Delta-vp.z),1e-6)
                self.assertLess(abs((sp.vx-p.vx)/Delta-vp.vx),1e-6)
                self.assertLess(abs((sp.vy-p.vy)/Delta-vp.vy),1e-6)
                self.assertLess(abs((sp.vz-p.vz)/Delta-vp.vz),1e-6)
                self.assertLess(abs((sp.m-p.m)/Delta-vp.m),1e-6)
                
    def test_all_2nd_order_init(self):
        vlist = ["a","e","i","Omega","omega","f","m"]
        # Testing a few random initial conditions
        paramslist = [ 
                (1e-3,1.,0.1,0.02,0.3,0.56,0.4),
                (1e-6,2.,0.02,0.0132,0.33,1.56,0.14),
                (234.3e-6,1.7567,0.561,0.572,0.573,2.56,0.354),
                (1e-2,1.7567,0.1561,0.15472,0.24573,12.56,1.354),
                (1e-7,3.7567,0.00061,0.23572,0.523473,2.56,3.354),
                ]
        for params in paramslist:
            for v1 in vlist:
                for v2 in vlist:
                    sim = rebound.Simulation()
                    sim.add(m=1.)
                    m,a,e,inc,Omega,omega,f= params
                    sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    sim.add(a=1.76)
                    p = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    vp = rebound.Particle(simulation=sim, primary=sim.particles[0],variation=v1,variation2=v2,variation_order=2,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    Delta = 1e-4
                    
                    if v1==v2:
                        m,a,e,inc,Omega,omega,f= params
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
                        m,a,e,inc,Omega,omega,f= params
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
                        prec = 1e-6
                        
                        self.assertLess(abs((sp.x-2.*p.x+sm.x)/Delta/Delta-vp.x),prec)
                        self.assertLess(abs((sp.y-2.*p.y+sm.y)/Delta/Delta-vp.y),prec)
                        self.assertLess(abs((sp.z-2.*p.z+sm.z)/Delta/Delta-vp.z),prec)
                        self.assertLess(abs((sp.vx-2.*p.vx+sm.vx)/Delta/Delta-vp.vx),prec)
                        self.assertLess(abs((sp.vy-2.*p.vy+sm.vy)/Delta/Delta-vp.vy),prec)
                        self.assertLess(abs((sp.vz-2.*p.vz+sm.vz)/Delta/Delta-vp.vz),prec)
                        self.assertLess(abs((sp.m-2.*p.m+sm.m)/Delta/Delta-vp.m),prec)
                    else:
                        m,a,e,inc,Omega,omega,f= params
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
                        
                        m,a,e,inc,Omega,omega,f= params
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

                        m,a,e,inc,Omega,omega,f= params
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

                        m,a,e,inc,Omega,omega,f= params
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

                        prec = 1e-6
                        self.assertLess(abs((spp.x  -spm.x  -smp.x  +smm.x  )/Delta/Delta/4.-vp.x),prec)
                        self.assertLess(abs((spp.y  -spm.y  -smp.y  +smm.y  )/Delta/Delta/4.-vp.y),prec)
                        self.assertLess(abs((spp.z  -spm.z  -smp.z  +smm.z  )/Delta/Delta/4.-vp.z),prec)
                        self.assertLess(abs((spp.vx -spm.vx -smp.vx +smm.vx )/Delta/Delta/4.-vp.vx),prec)
                        self.assertLess(abs((spp.vy -spm.vy -smp.vy +smm.vy )/Delta/Delta/4.-vp.vy),prec)
                        self.assertLess(abs((spp.vz -spm.vz -smp.vz +smm.vz )/Delta/Delta/4.-vp.vz),prec)
                        self.assertLess(abs((spp.m  -spm.m  -smp.m  +smm.m  )/Delta/Delta/4.-vp.m),prec)


class TestVariationalFull(unittest.TestCase):
    def test_all_1st_order_full(self):
        vlist = ["a","e","i","Omega","omega","f","m"]
        paramslist = [ 
                (1e-3,1.,0.1,0.02,0.3,0.56,0.4),
                (1e-6,2.,0.02,0.0132,0.33,1.56,0.14),
                (234.3e-6,1.7567,0.561,0.572,0.573,2.56,0.354),
                (1e-2,1.7567,0.1561,0.15472,0.24573,12.56,1.354),
                (1e-7,3.7567,0.00061,0.23572,0.523473,2.56,3.354),
                ]
        for params in paramslist:
            for v in vlist:
                m,a,e,inc,Omega,omega,f= params
                simvp = rebound.Simulation()
                simvp.add(m=1.)
                p = rebound.Particle(simulation=simvp, primary=simvp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                simvp.add(p)
                simvp.add(primary=simvp.particles[0],a=1.76)
                var_i = simvp.add_variation()
                vp = rebound.Particle(simulation=simvp, primary=simvp.particles[0],variation=v,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                var_i.particles[1] = vp
                simvp.integrate(10.)


                simsp = rebound.Simulation()
                simsp.add(m=1.)
                p = rebound.Particle(simulation=simsp, primary=simsp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                simsp.add(p)
                simsp.add(primary=simsp.particles[0],a=1.76)
                Delta=1e-8
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
                sp = rebound.Particle(simulation=simsp, primary=simsp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                simsp._particles[1] = sp
                simsp.integrate(10.)

                vp = var_i.particles[1]
                sp = simsp.particles[1]
                p = simvp.particles[1]

                prec = 1e-5
                self.assertLess(abs((sp.x-p.x)/Delta-vp.x),prec)
                self.assertLess(abs((sp.y-p.y)/Delta-vp.y),prec)
                self.assertLess(abs((sp.z-p.z)/Delta-vp.z),prec)
                self.assertLess(abs((sp.vx-p.vx)/Delta-vp.vx),prec)
                self.assertLess(abs((sp.vy-p.vy)/Delta-vp.vy),prec)
                self.assertLess(abs((sp.vz-p.vz)/Delta-vp.vz),prec)
                self.assertLess(abs((sp.m-p.m)/Delta-vp.m),prec)
    
    
    
    
    def test_all_2nd_order_full(self):
        self.run_2nd_order_full()
    def test_all_2nd_order_full_com(self):
        self.run_2nd_order_full(com=True)
    def run_2nd_order_full(self, com=False):
        vlist = ["a","e","i","Omega","omega","f","m"]
        # Testing a few random initial conditions
        paramslist = [ 
                (1e-3,1.,0.1,0.02,0.3,0.56,0.4),
                (1e-6,2.,0.02,0.0132,0.33,1.56,0.14),
                (234.3e-6,1.7567,0.561,0.572,0.573,2.56,0.354),
                (1e-2,1.7567,0.1561,0.15472,0.24573,12.56,1.354),
                (1e-7,3.7567,0.00061,0.23572,0.523473,2.56,3.354),
                ]
        for params in paramslist:
            for v1 in vlist:
                for v2 in vlist:
                    m,a,e,inc,Omega,omega,f= params
                    simvp = rebound.Simulation()
                    simvp.add(m=1.)
                    p = rebound.Particle(simulation=simvp, primary=simvp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    simvp.add(p)
                    simvp.add(primary=simvp.particles[0],a=1.76)
                    var_ia = simvp.add_variation()
                    var_ib = simvp.add_variation()
                    var_ii = simvp.add_variation(order=2,first_order=var_ia,first_order_2=var_ib)
                    var_ia.particles[1] = rebound.Particle(simulation=simvp, primary=simvp.particles[0],variation=v1,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    var_ib.particles[1] = rebound.Particle(simulation=simvp, primary=simvp.particles[0],variation=v2,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    vpii = rebound.Particle(simulation=simvp, primary=simvp.particles[0],variation=v1,variation2=v2,variation_order=2,m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                    var_ii.particles[1] = vpii
                    if com:
                        simvp.move_to_com()
                    simvp.integrate(10.)
                    
                    p = simvp.particles[1]
                    vp_ia = var_ia.particles[1]
                    vp_ib = var_ib.particles[1]
                    vp_ii = var_ii.particles[1]


                    Delta = 1e-4
                    
                    if v1==v2:
                        m,a,e,inc,Omega,omega,f= params
                        sim = rebound.Simulation()
                        sim.add(m=1.)
                        sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        sim.add(a=1.76)
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
                        sim._particles[1] = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            sim.move_to_com()
                        sim.integrate(10.)
                        sp = sim._particles[1]
                        
                        
                        m,a,e,inc,Omega,omega,f= params
                        sim = rebound.Simulation()
                        sim.add(m=1.)
                        sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        sim.add(a=1.76)
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
                        sim._particles[1] = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            sim.move_to_com()
                        sim.integrate(10.)
                        sm = sim._particles[1]
                        prec = 1e-4
                        
                        self.assertLess(abs((sp.x-2.*p.x+sm.x)/Delta/Delta-vp_ii.x),prec)
                        self.assertLess(abs((sp.y-2.*p.y+sm.y)/Delta/Delta-vp_ii.y),prec)
                        self.assertLess(abs((sp.z-2.*p.z+sm.z)/Delta/Delta-vp_ii.z),prec)
                        self.assertLess(abs((sp.vx-2.*p.vx+sm.vx)/Delta/Delta-vp_ii.vx),prec)
                        self.assertLess(abs((sp.vy-2.*p.vy+sm.vy)/Delta/Delta-vp_ii.vy),prec)
                        self.assertLess(abs((sp.vz-2.*p.vz+sm.vz)/Delta/Delta-vp_ii.vz),prec)
                        self.assertLess(abs((sp.m-2.*p.m+sm.m)/Delta/Delta-vp_ii.m),prec)
                        
                        # Also checking first order
                        prec = 1e-5
                        self.assertLess(abs((sp.x -sm.x )/Delta/2.-vp_ia.x ),prec)
                        self.assertLess(abs((sp.y -sm.y )/Delta/2.-vp_ia.y ),prec)
                        self.assertLess(abs((sp.z -sm.z )/Delta/2.-vp_ia.z ),prec)
                        self.assertLess(abs((sp.vx-sm.vx)/Delta/2.-vp_ia.vx),prec)
                        self.assertLess(abs((sp.vy-sm.vy)/Delta/2.-vp_ia.vy),prec)
                        self.assertLess(abs((sp.vz-sm.vz)/Delta/2.-vp_ia.vz),prec)
                        self.assertLess(abs((sp.m -sm.m) /Delta/2.-vp_ia.m ),prec)
                    else:
                        m,a,e,inc,Omega,omega,f= params
                        sim = rebound.Simulation()
                        sim.add(m=1.)
                        sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        sim.add(a=1.76)
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
                        sim._particles[1] = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            sim.move_to_com()
                        sim.integrate(10.)
                        spp = sim._particles[1]
                        
                        m,a,e,inc,Omega,omega,f= params
                        sim = rebound.Simulation()
                        sim.add(m=1.)
                        sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        sim.add(a=1.76)
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
                        sim._particles[1] = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            sim.move_to_com()
                        sim.integrate(10.)
                        spm = sim._particles[1]

                        m,a,e,inc,Omega,omega,f= params
                        sim = rebound.Simulation()
                        sim.add(m=1.)
                        sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        sim.add(a=1.76)
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
                        sim._particles[1] = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            sim.move_to_com()
                        sim.integrate(10.)
                        smp = sim._particles[1]

                        m,a,e,inc,Omega,omega,f= params
                        sim = rebound.Simulation()
                        sim.add(m=1.)
                        sim.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        sim.add(a=1.76)
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
                        sim._particles[1] = rebound.Particle(simulation=sim, primary=sim.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            sim.move_to_com()
                        sim.integrate(10.)
                        smm = sim._particles[1]

                        prec = 1e-4
                        self.assertLess(abs((spp.x  -spm.x  -smp.x  +smm.x  )/Delta/Delta/4.-vp_ii.x),prec)
                        self.assertLess(abs((spp.y  -spm.y  -smp.y  +smm.y  )/Delta/Delta/4.-vp_ii.y),prec)
                        self.assertLess(abs((spp.z  -spm.z  -smp.z  +smm.z  )/Delta/Delta/4.-vp_ii.z),prec)
                        self.assertLess(abs((spp.vx -spm.vx -smp.vx +smm.vx )/Delta/Delta/4.-vp_ii.vx),prec)
                        self.assertLess(abs((spp.vy -spm.vy -smp.vy +smm.vy )/Delta/Delta/4.-vp_ii.vy),prec)
                        self.assertLess(abs((spp.vz -spm.vz -smp.vz +smm.vz )/Delta/Delta/4.-vp_ii.vz),prec)
                        self.assertLess(abs((spp.m  -spm.m  -smp.m  +smm.m  )/Delta/Delta/4.-vp_ii.m),prec)

