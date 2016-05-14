import rebound
import unittest
import datetime


class TestVariational(unittest.TestCase):
    paramlist = [ 
            (1e-3,1.,0.1,0.02,0.3,0.56,0.4),
            (1e-6,2.,0.02,0.0132,0.33,1.56,0.14),
            (234.3e-6,1.7567,0.561,0.572,0.573,2.56,0.354),
            (1e-2,1.7567,0.1561,0.15472,0.24573,12.56,1.354),
            (1e-7,3.7567,0.00061,0.23572,0.523473,2.56,3.354),
            ]
    paramkeys = ["m","a","e","inc","Omega","omega","f"]
    def test_all_1st_order_full(self):
        for params in self.paramlist:
            for v in self.paramkeys:
                Delta=1e-8

                param = dict(zip(self.paramkeys, params))
                simvp = rebound.Simulation()
                simvp.add(m=1.)
                simvp.add(**param)
                simvp.add(primary=simvp.particles[0],a=1.76, m=1e-3)
                var_i = simvp.add_variation()
                var_i.vary(1,v)
                simvp.integrate(1.4)

                simsp = rebound.Simulation()
                simsp.add(m=1.)
                param[v] += Delta
                simsp.add(**param)
                simsp.add(primary=simsp.particles[0],a=1.76, m=1e-3)
                simsp.integrate(1.4)

                prec = 1e-5
                dp = (simsp.particles[1]-simvp.particles[1])/Delta - var_i.particles[1]
                self.assertLess(abs(dp.x ),prec)
                self.assertLess(abs(dp.y ),prec)
                self.assertLess(abs(dp.z ),prec)
                self.assertLess(abs(dp.vx),prec)
                self.assertLess(abs(dp.vy),prec)
                self.assertLess(abs(dp.vz),prec)
                self.assertLess(abs(dp.m ),prec)
    
    
    
    def test_all_2nd_order_full(self):
        self.run_2nd_order_full(com=False)
    def test_all_2nd_order_full_com(self):
        self.run_2nd_order_full(com=True)
    def run_2nd_order_full(self, com):
        for params in self.paramlist:
            for v1 in self.paramkeys:
                for v2 in self.paramkeys:
                    param = dict(zip(self.paramkeys, params))
                    simvp = rebound.Simulation()
                    simvp.add(m=1.)
                    simvp.add(**param)
                    simvp.add(primary=simvp.particles[0],a=1.76, m=1e-3)
                    var_ia = simvp.add_variation()
                    var_ib = simvp.add_variation()
                    var_ii = simvp.add_variation(order=2,first_order=var_ia,first_order_2=var_ib)
                    var_ia.vary(1, v1)
                    var_ib.vary(1, v2)
                    var_ii.vary(1, v1, v2)
                    if com:
                        simvp.move_to_com()
                    simvp.integrate(1.4)

                    Delta = 1e-5
                        
                    param = dict(zip(self.paramkeys, params))
                    simpp = rebound.Simulation()
                    simpp.add(m=1.)
                    param[v1] += Delta
                    param[v2] += Delta
                    simpp.add(**param)
                    simpp.add(primary=simpp.particles[0],a=1.76, m=1e-3)
                    if com:
                        simpp.move_to_com()
                    simpp.integrate(1.4)
                    
                    param = dict(zip(self.paramkeys, params))
                    simpm = rebound.Simulation()
                    simpm.add(m=1.)
                    param[v1] += Delta
                    param[v2] -= Delta
                    simpm.add(**param)
                    simpm.add(primary=simpm.particles[0],a=1.76, m=1e-3)
                    if com:
                        simpm.move_to_com()
                    simpm.integrate(1.4)

                    param = dict(zip(self.paramkeys, params))
                    simmp = rebound.Simulation()
                    simmp.add(m=1.)
                    param[v1] -= Delta
                    param[v2] += Delta
                    simmp.add(**param)
                    simmp.add(primary=simmp.particles[0],a=1.76, m=1e-3)
                    if com:
                        simmp.move_to_com()
                    simmp.integrate(1.4)

                    param = dict(zip(self.paramkeys, params))
                    simmm = rebound.Simulation()
                    simmm.add(m=1.)
                    param[v1] -= Delta
                    param[v2] -= Delta
                    simmm.add(**param)
                    simmm.add(primary=simmm.particles[0],a=1.76, m=1e-3)
                    if com:
                        simmm.move_to_com()
                    simmm.integrate(1.4)

                    prec = 1e-4
                    
                    dp = (simpp.particles[1]-simpm.particles[1]-simmp.particles[1]+simmm.particles[1])/(Delta*Delta*4.) - var_ii.particles[1]

                    self.assertLess(abs(dp.x),prec)
                    self.assertLess(abs(dp.y),prec)
                    self.assertLess(abs(dp.z),prec)
                    self.assertLess(abs(dp.vx),prec)
                    self.assertLess(abs(dp.vy),prec)
                    self.assertLess(abs(dp.vz),prec)
                    self.assertLess(abs(dp.m),prec)

class TestVariationalPal(TestVariational):
    paramkeys = ["m","a","h","k","l","ix","iy"]
    paramlist = [ 
            [1e-3,      1.,     0.,     0.,     0.,         0.0,     0.0],
            [1e-3,      1.,     0.1,    0.02,   0.3,        0.0,     0.0],
            [1e-6,      2.,     0.02,   0.0132, 0.33,       0.126,   0.14],
            [234.3e-6,  1.7567, 0.561,  0.572,  0.573,      0.056,   0.0091354],
            [1e-2,      1.7567, 0.1561, 0.15472,0.24573,    0.0056,  0.0013],
            [1e-7,      3.7567, 0.00061,0.23572,0.523473,   0.47256, 0.000024],
            [1e-7,      3.7567, 0.00061,0.23572,0.523473,   1.97,    0.0],
            [1e-7,      3.7567, 0.00061,0.23572,0.523473,   0.0,     1.97],
            ]

    
class TestVariationalTestParticle(unittest.TestCase):
    paramlist = [ 
            (1.,0.1,0.02,0.3,0.56,0.4),
            (2.,0.02,0.0132,0.33,1.56,0.14),
            (1.7567,0.561,0.572,0.573,2.56,0.354),
            (1.7567,0.1561,0.15472,0.24573,12.56,1.354),
            (3.7567,0.00061,0.23572,0.523473,2.56,3.354),
            ]
    paramkeys = ["a","e","inc","Omega","omega","f"]
    def test_all_1st_order_full(self):
        for params in self.paramlist:
            for v in self.paramkeys:
                Delta=1e-8

                param = dict(zip(self.paramkeys, params))
                simvp = rebound.Simulation()
                simvp.add(m=1.)
                simvp.add(**param)
                simvp.add(primary=simvp.particles[0],a=1.76, m=1e-3)
                var_i = simvp.add_variation(testparticle=1)
                var_i.vary(1,v)
                simvp.integrate(1.4)

                simsp = rebound.Simulation()
                simsp.add(m=1.)
                param[v] += Delta
                simsp.add(**param)
                simsp.add(primary=simsp.particles[0],a=1.76, m=1e-3)
                simsp.integrate(1.4)

                prec = 1e-5
                dp = (simsp.particles[1]-simvp.particles[1])/Delta - var_i.particles[0]
                self.assertLess(abs(dp.x ),prec)
                self.assertLess(abs(dp.y ),prec)
                self.assertLess(abs(dp.z ),prec)
                self.assertLess(abs(dp.vx),prec)
                self.assertLess(abs(dp.vy),prec)
                self.assertLess(abs(dp.vz),prec)
                self.assertLess(abs(dp.m ),prec)
    
    
    
    def test_all_2nd_order_full(self):
        self.run_2nd_order_full(com=False)
    def test_all_2nd_order_full_com(self):
        self.run_2nd_order_full(com=True)
    def run_2nd_order_full(self, com):
        for params in self.paramlist:
            for v1 in self.paramkeys:
                for v2 in self.paramkeys:
                    param = dict(zip(self.paramkeys, params))
                    simvp = rebound.Simulation()
                    simvp.add(m=1.)
                    simvp.add(**param)
                    simvp.add(primary=simvp.particles[0],a=1.76, m=1e-3)
                    var_ia = simvp.add_variation(testparticle=1)
                    var_ib = simvp.add_variation(testparticle=1)
                    var_ii = simvp.add_variation(order=2,first_order=var_ia,first_order_2=var_ib,testparticle=1)
                    var_ia.vary(1, v1)
                    var_ib.vary(1, v2)
                    var_ii.vary(1, v1, v2)
                    if com:
                        simvp.move_to_com()
                    simvp.integrate(1.4)

                    Delta = 1e-5
                        
                    param = dict(zip(self.paramkeys, params))
                    simpp = rebound.Simulation()
                    simpp.add(m=1.)
                    param[v1] += Delta
                    param[v2] += Delta
                    simpp.add(**param)
                    simpp.add(primary=simpp.particles[0],a=1.76, m=1e-3)
                    if com:
                        simpp.move_to_com()
                    simpp.integrate(1.4)
                    
                    param = dict(zip(self.paramkeys, params))
                    simpm = rebound.Simulation()
                    simpm.add(m=1.)
                    param[v1] += Delta
                    param[v2] -= Delta
                    simpm.add(**param)
                    simpm.add(primary=simpm.particles[0],a=1.76, m=1e-3)
                    if com:
                        simpm.move_to_com()
                    simpm.integrate(1.4)

                    param = dict(zip(self.paramkeys, params))
                    simmp = rebound.Simulation()
                    simmp.add(m=1.)
                    param[v1] -= Delta
                    param[v2] += Delta
                    simmp.add(**param)
                    simmp.add(primary=simmp.particles[0],a=1.76, m=1e-3)
                    if com:
                        simmp.move_to_com()
                    simmp.integrate(1.4)

                    param = dict(zip(self.paramkeys, params))
                    simmm = rebound.Simulation()
                    simmm.add(m=1.)
                    param[v1] -= Delta
                    param[v2] -= Delta
                    simmm.add(**param)
                    simmm.add(primary=simmm.particles[0],a=1.76, m=1e-3)
                    if com:
                        simmm.move_to_com()
                    simmm.integrate(1.4)

                    prec = 1e-4
                    
                    dp = (simpp.particles[1]-simpm.particles[1]-simmp.particles[1]+simmm.particles[1])/(Delta*Delta*4.) - var_ii.particles[0]

                    self.assertLess(abs(dp.x),prec)
                    self.assertLess(abs(dp.y),prec)
                    self.assertLess(abs(dp.z),prec)
                    self.assertLess(abs(dp.vx),prec)
                    self.assertLess(abs(dp.vy),prec)
                    self.assertLess(abs(dp.vz),prec)
                    self.assertLess(abs(dp.m),prec)

class TestVariationalPalTestParticle(TestVariationalTestParticle):
    paramkeys = ["a","h","k","l","ix","iy"]
    paramlist = [ 
            [1.,     0.,     0.,     0.,         0.0,     0.0],
            [1.,     0.1,    0.02,   0.3,        0.0,     0.0],
            [2.,     0.02,   0.0132, 0.33,       0.126,   0.14],
            [1.7567, 0.561,  0.572,  0.573,      0.056,   0.0091354],
            [1.7567, 0.1561, 0.15472,0.24573,    0.0056,  0.0013],
            [3.7567, 0.00061,0.23572,0.523473,   0.47256, 0.000024],
            [3.7567, 0.00061,0.23572,0.523473,   1.97,    0.0],
            [3.7567, 0.00061,0.23572,0.523473,   0.0,     1.97],
            ]

if __name__ == "__main__":
    unittest.main()
