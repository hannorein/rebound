import rebound
import unittest
import datetime

vlist = ["m","a","h","k","l","ix","iy"]
paramslist = [ 
        [1e-3,      1.,     0.,     0.,     0.,         0.0,     0.0],
        [1e-3,      1.,     0.1,    0.02,   0.3,        0.0,     0.0],
        [1e-6,      2.,     0.02,   0.0132, 0.33,       0.126,   0.14],
        [234.3e-6,  1.7567, 0.561,  0.572,  0.573,      0.056,   0.0091354],
        [1e-2,      1.7567, 0.1561, 0.15472,0.24573,    0.0056,  0.0013],
        [1e-7,      3.7567, 0.00061,0.23572,0.523473,   0.47256, 0.000024],
        [1e-7,      3.7567, 0.00061,0.23572,0.523473,   1.97,    0.0],
        [1e-7,      3.7567, 0.00061,0.23572,0.523473,   0.0,     1.97],
        ]

class TestVariationalPal(unittest.TestCase):

    def get_sim(self, params):
        args = dict(zip(vlist,params))
        sim = rebound.Simulation()
        sim.add(m=1.)
        p = rebound.Particle(simulation=sim, primary=sim.particles[0], **args)
        sim.add(p)
        sim.add(primary=sim.particles[0],a=1.76)
        sim.integrate(10.)
        return sim


    def test_all_1st_order(self):
        Delta = 1e-8
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
                simsp = self.get_sim(params)
                params[vi] -= Delta

                vp = var_i.particles[1]
                sp = simsp.particles[1]
                p = simvp.particles[1]

                p_fd = (sp-p)/Delta

                prec = 1e-5
                msg = "Problem with '%s' derivative for the following parameters: %s." % (v,args)
                self.assertLess(abs(p_fd.x -vp.x ), prec, msg=msg)
                self.assertLess(abs(p_fd.y -vp.y ), prec, msg=msg)
                self.assertLess(abs(p_fd.z -vp.z ), prec, msg=msg)
                self.assertLess(abs(p_fd.vx-vp.vx), prec, msg=msg)
                self.assertLess(abs(p_fd.vy-vp.vy), prec, msg=msg)
                self.assertLess(abs(p_fd.vz-vp.vz), prec, msg=msg)
                self.assertLess(abs(p_fd.m -vp.m ), prec, msg=msg)
    
    def test_all_2nd_order(self):
        Delta = 1e-5
        for params in paramslist:
            for vi,v in enumerate(vlist):
                for wi,w in enumerate(vlist):
                    if (w=="m" or v=="m") and params[0]<1e-4:
                        # Ignore mass derivates for very small particles (finite differences don't work).
                        continue
                    args = dict(zip(vlist,params))
                    simvp = rebound.Simulation()
                    simvp.add(m=1.)
                    p = rebound.Particle(simulation=simvp, primary=simvp.particles[0], **args)
                    simvp.add(p)
                    simvp.add(primary=simvp.particles[0],a=1.76)
                    var_v = simvp.add_variation()
                    var_v.vary_pal(1,v)
                    var_w = simvp.add_variation()
                    var_w.vary_pal(1,w)
                    var_vw = simvp.add_variation(order=2, first_order=var_w, first_order_2=var_v)
                    var_vw.vary_pal(1,w,v)
                    simvp.integrate(10.)

                    params[vi] += Delta
                    params[wi] += Delta
                    simspp = self.get_sim(params)
                    params[vi] -= Delta
                    params[wi] -= Delta

                    params[vi] += Delta
                    params[wi] -= Delta
                    simspm = self.get_sim(params)
                    params[vi] -= Delta
                    params[wi] += Delta
                    
                    params[vi] -= Delta
                    params[wi] += Delta
                    simsmp = self.get_sim(params)
                    params[vi] += Delta
                    params[wi] -= Delta
                    
                    params[vi] -= Delta
                    params[wi] -= Delta
                    simsmm = self.get_sim(params)
                    params[vi] += Delta
                    params[wi] += Delta

                    vp  = var_vw.particles[1]
                    spp = simspp.particles[1]
                    spm = simspm.particles[1]
                    smp = simsmp.particles[1]
                    smm = simsmm.particles[1]

                    p_fd = (spp-spm-smp+smm)/(4.*Delta*Delta)

                    prec = 1e-4
                    msg ="Problem with '%s_%s' derivative for the following parameters: %s." % (v,w,args)
                    self.assertLess(abs(p_fd.x -vp.x ), prec, msg=msg)
                    self.assertLess(abs(p_fd.y -vp.y ), prec, msg=msg)
                    self.assertLess(abs(p_fd.z -vp.z ), prec, msg=msg)
                    self.assertLess(abs(p_fd.vx-vp.vx), prec, msg=msg)
                    self.assertLess(abs(p_fd.vy-vp.vy), prec, msg=msg)
                    self.assertLess(abs(p_fd.vz-vp.vz), prec, msg=msg)
                    self.assertLess(abs(p_fd.m -vp.m ), prec, msg=msg)

if __name__ == "__main__":
    unittest.main()
