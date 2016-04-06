import rebound
import unittest
import datetime


class TestVariationalPal(unittest.TestCase):
    def test_all_1st_order_full(self):
        Delta = 1e-8
        vlist = ["m","a","h","k","l","ix","iy"]
        paramslist = [ 
                [1e-3,      1.,     0.1,    0.02,   0.3,        0.0,     0.0],
                [1e-6,      2.,     0.02,   0.0132, 0.33,       0.126,   0.14],
                [234.3e-6,  1.7567, 0.561,  0.572,  0.573,      0.056,   0.0091354],
                [1e-2,      1.7567, 0.1561, 0.15472,0.24573,    0.0056,  0.0013],
                [1e-7,      3.7567, 0.00061,0.23572,0.523473,   0.47256, 0.000024],
                [1e-7,      3.7567, 0.00061,0.23572,0.523473,   2., 0.0],
                ]
        for params in paramslist:
            for vi,v in enumerate(vlist):
                args = dict(zip(vlist,params))
                simvp = rebound.Simulation()
                simvp.add(m=1.)
                p = rebound.Particle(simulation=simvp, primary=simvp.particles[0], **args)
                simvp.add(p)
                simvp.add(primary=simvp.particles[0],a=1.76)
                var_i = simvp.add_variation()
                var_i.vary_pal(1,v)
                simvp.integrate(10.)

                params[vi] += Delta
                args = dict(zip(vlist,params))
                simsp = rebound.Simulation()
                simsp.add(m=1.)
                p = rebound.Particle(simulation=simsp, primary=simsp.particles[0], **args)
                simsp.add(p)
                simsp.add(primary=simsp.particles[0],a=1.76)
                simsp.integrate(10.)
                params[vi] -= Delta

                vp = var_i.particles[1]
                sp = simsp.particles[1]
                p = simvp.particles[1]

                p_fd = (sp-p)/Delta

                prec = 1e-5
                self.assertLess(abs(p_fd.x -vp.x ), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
                self.assertLess(abs(p_fd.y -vp.y ), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
                self.assertLess(abs(p_fd.z -vp.z ), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
                self.assertLess(abs(p_fd.vx-vp.vx), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
                self.assertLess(abs(p_fd.vy-vp.vy), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
                self.assertLess(abs(p_fd.vz-vp.vz), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
                self.assertLess(abs(p_fd.m -vp.m ), prec, msg="Problem with '%s' derivative for the following parameters: %s."%(v,args))
    
    
    
if __name__ == "__main__":
    unittest.main()
