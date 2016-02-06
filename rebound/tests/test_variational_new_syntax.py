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
                var_i = simvp.add_variation()
                var_i.vary(1,v)
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
                simsp.particles[1] = sp
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
                    var_ia.vary(1, v1)
                    var_ib.vary(1, v2)
                    var_ii.vary(1, v1, v2)
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
                        simp = rebound.Simulation()
                        simp.add(m=1.)
                        simp.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        simp.add(a=1.76)
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
                        simp.particles[1] = rebound.Particle(simulation=simp, primary=simp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            simp.move_to_com()
                        simp.integrate(10.)
                        sp = simp.particles[1]
                        
                        
                        m,a,e,inc,Omega,omega,f= params
                        simm = rebound.Simulation()
                        simm.add(m=1.)
                        simm.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        simm.add(a=1.76)
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
                        simm.particles[1] = rebound.Particle(simulation=simm, primary=simm.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            simm.move_to_com()
                        simm.integrate(10.)
                        sm = simm.particles[1]
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
                        simpp = rebound.Simulation()
                        simpp.add(m=1.)
                        simpp.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        simpp.add(a=1.76)
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
                        simpp.particles[1] = rebound.Particle(simulation=simpp, primary=simpp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            simpp.move_to_com()
                        simpp.integrate(10.)
                        spp = simpp.particles[1]
                        
                        m,a,e,inc,Omega,omega,f= params
                        simpm = rebound.Simulation()
                        simpm.add(m=1.)
                        simpm.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        simpm.add(a=1.76)
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
                        simpm.particles[1] = rebound.Particle(simulation=simpm, primary=simpm.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            simpm.move_to_com()
                        simpm.integrate(10.)
                        spm = simpm.particles[1]

                        m,a,e,inc,Omega,omega,f= params
                        simmp = rebound.Simulation()
                        simmp.add(m=1.)
                        simmp.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        simmp.add(a=1.76)
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
                        simmp.particles[1] = rebound.Particle(simulation=simmp, primary=simmp.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            simmp.move_to_com()
                        simmp.integrate(10.)
                        smp = simmp.particles[1]

                        m,a,e,inc,Omega,omega,f= params
                        simmm = rebound.Simulation()
                        simmm.add(m=1.)
                        simmm.add(m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        simmm.add(a=1.76)
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
                        simmm.particles[1] = rebound.Particle(simulation=simmm, primary=simmm.particles[0],m=m,a=a,e=e,inc=inc,Omega=Omega,omega=omega,f=f)
                        if com:
                            simmm.move_to_com()
                        simmm.integrate(10.)
                        smm = simmm.particles[1]

                        prec = 1e-4
                        self.assertLess(abs((spp.x  -spm.x  -smp.x  +smm.x  )/Delta/Delta/4.-vp_ii.x),prec)
                        self.assertLess(abs((spp.y  -spm.y  -smp.y  +smm.y  )/Delta/Delta/4.-vp_ii.y),prec)
                        self.assertLess(abs((spp.z  -spm.z  -smp.z  +smm.z  )/Delta/Delta/4.-vp_ii.z),prec)
                        self.assertLess(abs((spp.vx -spm.vx -smp.vx +smm.vx )/Delta/Delta/4.-vp_ii.vx),prec)
                        self.assertLess(abs((spp.vy -spm.vy -smp.vy +smm.vy )/Delta/Delta/4.-vp_ii.vy),prec)
                        self.assertLess(abs((spp.vz -spm.vz -smp.vz +smm.vz )/Delta/Delta/4.-vp_ii.vz),prec)
                        self.assertLess(abs((spp.m  -spm.m  -smp.m  +smm.m  )/Delta/Delta/4.-vp_ii.m),prec)

