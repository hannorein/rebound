import rebound
import unittest
import datetime


class TestVariationalFull_NewSyntax(unittest.TestCase):
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
                var_i = simvp.add_variational()
                simvp.init_variational_particle(var_i,1,v)
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

                vp = simvp.particles[var_i+1]
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
                    var_ia = simvp.add_variational()
                    var_ib = simvp.add_variational()
                    var_ii = simvp.add_variational(order=2,index_1st_order_a=var_ia,index_1st_order_b=var_ib)
                    simvp.init_variational_particle(var_ia, 1, v1)
                    simvp.init_variational_particle(var_ib, 1, v2)
                    simvp.init_variational_particle(var_ii, 1, v1, v2)
                    if com:
                        simvp.move_to_com()
                    simvp.integrate(10.)
                    
                    p = simvp.particles[1]
                    vp_ia = simvp.particles[var_ia+1]
                    vp_ib = simvp.particles[var_ib+1]
                    vp_ii = simvp.particles[var_ii+1]


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

