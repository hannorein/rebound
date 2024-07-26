import rebound
import unittest
import math
import rebound.data


def derivatives_ho(ode, yDot, y, t):
    m = 1.
    k = 100.
    yDot[0] = y[1]
    yDot[1] = -k/m*y[0]

class TestIntegratorBSHarmonic(unittest.TestCase):
    def test_bs_harmonic_only(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        ode_ho = sim.create_ode(length=2, needs_nbody=False)
        ode_ho.derivatives = derivatives_ho

        ode_ho.y[0] = 1. 
        ode_ho.y[1] = 0. # zero velocity

        sim.integrate(20.*math.pi)
        self.assertLess(math.fabs(ode_ho.y[0]-1.),2e-10)
        self.assertLess(math.fabs(ode_ho.y[1]),2e-9)
   
    def test_bs_harmonic_only_low_eps(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.ri_bs.eps_abs = 1e-5
        sim.ri_bs.eps_rel = 1e-5
        ode_ho = sim.create_ode(length=2, needs_nbody=False)
        ode_ho.derivatives = derivatives_ho

        ode_ho.y[0] = 1. 
        ode_ho.y[1] = 0. # zero velocity

        sim.integrate(20.*math.pi)
        self.assertLess(math.fabs(ode_ho.y[0]-1.),2e-8)
        self.assertLess(math.fabs(ode_ho.y[1]),5e-6)
    
    def test_bs_harmonic_only_high_eps(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.ri_bs.eps_abs = 1e-10
        sim.ri_bs.eps_rel = 1e-10
        ode_ho = sim.create_ode(length=2, needs_nbody=False)
        ode_ho.derivatives = derivatives_ho

        ode_ho.y[0] = 1. 
        ode_ho.y[1] = 0. # zero velocity

        sim.integrate(20.*math.pi)
        self.assertLess(math.fabs(ode_ho.y[0]-1.),1e-10)
        self.assertLess(math.fabs(ode_ho.y[1]),1e-10)
    
    
    def test_bs_harmonic_with_nbody(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.123);
        sim.add(m=1e-3,a=2.6,e=0.123);
        sim.integrator = "BS"
        ode_ho = sim.create_ode(length=2, needs_nbody=False)
        ode_ho.derivatives = derivatives_ho

        ode_ho.y[0] = 1. 
        ode_ho.y[1] = 0. # zero velocity

        sim.integrate(20.*math.pi)
        self.assertLess(math.fabs(ode_ho.y[0]-1.),2e-10)
        self.assertLess(math.fabs(ode_ho.y[1]),2e-9)
    
    def test_bs_harmonic_with_nbody_coupledy(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.123);
        sim.add(m=1e-3,a=2.6,e=0.123);
        sim.integrator = "BS"
        ode_ho = sim.create_ode(length=2, needs_nbody=True)
        ode_ho.derivatives = derivatives_ho

        ode_ho.y[0] = 1. 
        ode_ho.y[1] = 0. # zero velocity

        sim.integrate(20.*math.pi)
        self.assertLess(math.fabs(ode_ho.y[0]-1.),2e-10)
        self.assertLess(math.fabs(ode_ho.y[1]),2e-9)


    
def af(simp):
    sim = simp.contents
    x = sim.particles[0].x
    y = sim.particles[0].y
    z = sim.particles[0].z
    r = math.sqrt(x*x+y*y+z*z)
    sim.particles[0].ax -= x/(r*r*r)
    sim.particles[0].ay -= y/(r*r*r)
    sim.particles[0].az -= z/(r*r*r)
    
class TestIntegratorBS(unittest.TestCase):
    def test_bs_additional_force_only(self):
        sim = rebound.Simulation()
        sim.additional_forces = af
        sim.integrator = "bs"
        eps = 1e-11
        sim.ri_bs.eps_rel = eps
        sim.ri_bs.eps_abs = eps
        sim.add(m=0,x=1,vy=1)
        sim.integrate(2.*math.pi)
        self.assertLess(math.fabs(sim.particles[0].x-1.),5*eps)
        self.assertLess(math.fabs(sim.particles[0].vy-1.),5*eps)
        self.assertLess(math.fabs(sim.particles[0].y),5*eps)
        self.assertLess(math.fabs(sim.particles[0].vx),5*eps)


    def test_bs_outersolarsystem(self):
        for eps in [1e-5, 1e-7, 1e-9, 1e-11]:
            sim = rebound.Simulation()
            rebound.data.add_outer_solar_system(sim)
            sim.integrator = "bs"
            sim.ri_bs.eps_rel = eps
            sim.ri_bs.eps_abs = eps
            e0 = sim.energy()
            sim.integrate(1000)
            e1 = sim.energy()
            self.assertLess(math.fabs((e0-e1)/e1),5*eps)

    def test_bs_high_e(self):
        for eps in [1e-7, 1e-9, 1e-11]:
            sim = rebound.Simulation()
            sim.add(m=1)
            sim.add(m=1e-3,a=1,e=0.9)
            sim.add(m=1e-3,a=6,e=0.9,f=0.5,omega=1.6)
            sim.integrator = "bs"
            sim.ri_bs.eps_rel = eps
            sim.ri_bs.eps_abs = eps
            e0 = sim.energy()
            sim.integrate(1000)
            e1 = sim.energy()
            self.assertLess(math.fabs((e0-e1)/e1),60*eps)
    
    def test_bs_inout(self):
        sim = rebound.Simulation()
        rebound.data.add_outer_solar_system(sim)
        eps = 1e-6
        sim.integrator = "bs"
        sim.ri_bs.eps_rel = eps
        sim.ri_bs.eps_abs = eps
        sim.save_to_file("sim0.bin")
        sim1 = rebound.Simulation("sim0.bin")
        sim.integrate(100)
        sim1.integrate(100)
        sim1.save_to_file("sim1.bin")
        sim2 = rebound.Simulation("sim1.bin")
        sim.integrate(200)
        sim1.integrate(200)
        sim2.integrate(200)
        self.assertEqual(sim.particles[1].x, sim1.particles[1].x)
        self.assertEqual(sim.particles[1].x, sim2.particles[1].x)
        self.assertEqual(sim.particles[2].vx, sim1.particles[2].vx)
        self.assertEqual(sim.particles[2].vx, sim2.particles[2].vx)
        self.assertEqual(sim.t, sim1.t)
        self.assertEqual(sim.t, sim2.t)


    def test_bs_archive(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1)
        sim.add(m=1e-3,a=2,e=0.1)
        sim.save_to_file("test.sa",interval=10, delete_file=True)
        sim.integrate(100, exact_finish_time=0)
        sim1 = rebound.Simulationarchive("test.sa")[-1]
        sim.integrate(200, exact_finish_time=0)
        sim1.integrate(200, exact_finish_time=0)
        self.assertEqual(sim.particles[1].x, sim1.particles[1].x)
        self.assertEqual(sim.particles[2].vx, sim1.particles[2].vx)
        self.assertEqual(sim.t, sim1.t)
    
    def test_bs_change_particle_N(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.add(m=1)
        sim.add(m=1e-3,a=1,e=0.1)
        sim.N_active = sim.N
        for i in range(100):
            sim.add(a=1.4+i*0.01)
        sim.integrate(10)
        for i in range(100):
            sim.add(a=2.4+i*0.01)
        sim.integrate(20)
        for i in range(50):
            sim.remove(i*2+2)
        sim.integrate(20)
        self.assertEqual(sim.N, 150+2)
    
    def test_bs_collide(self):
        sim = rebound.Simulation()
        sim.integrator = "BS"
        sim.add(m=1,r=1,x=-3)
        sim.add(m=1,r=1,x=3)
        sim.collision = "direct"
        sim.collision_resolve = "merge"
        sim.integrate(20)
        self.assertEqual(sim.N, 1)
        self.assertEqual(sim.particles[0].m, 2.)
        self.assertEqual(sim.particles[0].x, 0.)
        self.assertEqual(sim.particles[0].vx, 0.)


class TestVariationalBS(unittest.TestCase):
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
                simvp.integrator = "BS"
                simvp.ri_bs.eps_abs = 1e-15
                simvp.ri_bs.eps_rel = 1e-15
                simvp.add(m=1.)
                simvp.add(**param)
                simvp.add(primary=simvp.particles[0],a=1.76, m=1e-3)
                var_i = simvp.add_variation()
                var_i.vary(1,v)
                simvp.integrate(1.4)

                simsp = rebound.Simulation()
                simsp.integrator = "BS"
                simsp.ri_bs.eps_abs = 1e-15
                simsp.ri_bs.eps_rel = 1e-15
                simsp.add(m=1.)
                param[v] += Delta
                simsp.add(**param)
                simsp.add(primary=simsp.particles[0],a=1.76, m=1e-3)
                simsp.integrate(1.4)

                prec = 2e-5
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
                    simvp.integrator = "BS"
                    simvp.ri_bs.eps_abs = 1e-15
                    simvp.ri_bs.eps_rel = 1e-15
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
                    simpp.integrator = "BS"
                    simpp.ri_bs.eps_abs = 1e-15
                    simpp.ri_bs.eps_rel = 1e-15
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
                    simpm.integrator = "BS"
                    simpm.ri_bs.eps_abs = 1e-15
                    simpm.ri_bs.eps_rel = 1e-15
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
                    simmp.integrator = "BS"
                    simmp.ri_bs.eps_abs = 1e-15
                    simmp.ri_bs.eps_rel = 1e-15
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
                    simmm.integrator = "BS"
                    simmm.ri_bs.eps_abs = 1e-15
                    simmm.ri_bs.eps_rel = 1e-15
                    simmm.add(m=1.)
                    param[v1] -= Delta
                    param[v2] -= Delta
                    simmm.add(**param)
                    simmm.add(primary=simmm.particles[0],a=1.76, m=1e-3)
                    if com:
                        simmm.move_to_com()
                    simmm.integrate(1.4)

                    prec = 5e-4
                    
                    dp = (simpp.particles[1]-simpm.particles[1]-simmp.particles[1]+simmm.particles[1])/(Delta*Delta*4.) - var_ii.particles[1]

                    self.assertLess(abs(dp.x),prec)
                    self.assertLess(abs(dp.y),prec)
                    self.assertLess(abs(dp.z),prec)
                    self.assertLess(abs(dp.vx),prec)
                    self.assertLess(abs(dp.vy),prec)
                    self.assertLess(abs(dp.vz),prec)
                    self.assertLess(abs(dp.m),prec)

class TestVariationalPal(TestVariationalBS):
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


if __name__ == "__main__":
    unittest.main()
