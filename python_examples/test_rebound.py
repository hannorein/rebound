import unittest
import rebound
import random
import numpy

class cartesian_to_orbital(unittest.TestCase):
    sun = rebound.Particle( m=1.,x=0.,y=0.,z=0.,vx=0.,vy=0.,vz=0.)
    G = 1.
    N_random_tests = 10 # num of random orbits to test in test_rand_r_to_orb_(f or M)
    
    def test_keplers_eq(self):
        '''test Kepler's equation'''
        zero_to_pi = numpy.linspace(0,numpy.pi,10,endpoint=True)
        pi_to_2pi = numpy.linspace(numpy.pi,2*numpy.pi,10,endpoint=False)
        e = numpy.asarray([0.,0.7,0.9,0.999])
        for ecc in e:
            # E & M match at 0,pi, so if M is in range 0<M<pi (or pi<M<2pi), so should E
            # Also check that Kepler's equation is satisfied to abs precision of 1e-15
            for M in zero_to_pi:
                E = rebound.get_E(ecc,M)
                err = E-ecc*numpy.sin(E)-M
                self.assertTrue(0.<= E <= numpy.pi, 
                                "E & M not in same half-plane: E=%.3f, M=%.3f"%(E,M))
                self.assertAlmostEqual(numpy.fabs(err),0.,places=14,msg=
                "Kepler's Eq not satisfied: e=%.3f, M=%.3f, abs. error=%.3e"%(ecc,M,err))
            for M in pi_to_2pi:
                E = rebound.get_E(ecc,M)
                err = E-ecc*numpy.sin(E)-M
                self.assertTrue(numpy.pi<= E <= 2*numpy.pi,
                                "E & M not in same half-plane: E=%.3f, M=%.3f"%(E,M))       
                self.assertAlmostEqual(numpy.fabs(err),0.,places=14,msg=
                "Kepler's Eq not satisfied: e=%.3f, M=%.3f, abs. error=%.3e"%(ecc,M,err))

    def test_r_to_orb_defaults(self):
        '''test conversion from orbital elements to cartesian and back 
        when not all orbital elements are passed to init_planet'''
        specified_cases = ((self.sun,1.,),
                       (self.sun,1.,0.01)
                       )
        results = ((1.,0.,0.,0.,0.,0.),
                   (1.,0.01,0.,0.,0.,0.))
        for j in range(len(specified_cases)):
            p = rebound.init_planet(*specified_cases[j])
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(results[j][0],o.a,places=12)
            self.assertAlmostEqual(results[j][1],o.e,places=12)
            self.assertAlmostEqual(results[j][2],o.inc,places=12)
            self.assertAlmostEqual(results[j][3],rebound.mod2pi(o.Omega),places=12)
            self.assertAlmostEqual(results[j][4],rebound.mod2pi(o.omega),places=12)
            self.assertAlmostEqual(results[j][5],o.f,places=12)
            
    def test_specified_r_to_orb(self):
        '''test conversion from orbital elements to cartesian and back with specified cases'''
        specified_cases = ((self.sun,1.,0.01,0.,0.,0.,0.,0.),
                       #(self.sun,1.,0.999,2.,3.,3.,0.,0.),
                       #(self.sun,1.759, 0.851, 0.882, 5.852, 6.139, 5.445,0.),
                       )
        for params in specified_cases:
            #print("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f"%(params[3:9]))
            p = rebound.init_planet(*params)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(params[1],o.a,places=12)
            self.assertAlmostEqual(params[2],o.e,places=12)
            self.assertAlmostEqual(params[5],o.inc,places=12)
            self.assertAlmostEqual(params[6],rebound.mod2pi(o.Omega),places=12)
            self.assertAlmostEqual(params[3],rebound.mod2pi(o.omega),places=12)
            self.assertAlmostEqual(params[4],o.f,places=12)   
            
    def test_rand_r_to_orb_f(self):    
        '''test conversion from orb. elements to cart. and back with random cases w/ true anom'''
        for q in range(self.N_random_tests):
            _a=random.uniform(1.,2.)
            _e=random.uniform(0.,1.)
            _inc=random.uniform(0,numpy.pi)
            _Omega=random.uniform(0,2*numpy.pi)
            _omega=random.uniform(0,2*numpy.pi)
            _f=random.uniform(0,2*numpy.pi)
            #print("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f"%(_a,_e,_inc,_Omega,_omega,_f))
            p = rebound.init_planet(self.sun,a=_a, e=_e, omega=_omega, f=_f, inc=_inc, 
                                    Omega=_Omega, mass=0.)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(_a,o.a,places=12)
            self.assertAlmostEqual(_e,o.e,places=12)
            self.assertAlmostEqual(_inc,o.inc,places=12)
            self.assertAlmostEqual(_Omega,rebound.mod2pi(o.Omega),places=12)
            self.assertAlmostEqual(_omega,rebound.mod2pi(o.omega),places=12)
            self.assertAlmostEqual(_f,rebound.mod2pi(o.f),places=12)
            
    def test_rand_r_to_orb_f(self):    
        '''test conversion from orb. elements to cart. and back with random cases w/ mean anom'''
        for q in range(self.N_random_tests):
            _a=random.uniform(1.,2.)
            _e=random.uniform(0.,1.)
            _inc=random.uniform(0,numpy.pi)
            _Omega=random.uniform(0,2*numpy.pi)
            _omega=random.uniform(0,2*numpy.pi)
            _M=random.uniform(0,2*numpy.pi)
            #print("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f"%(_a,_e,_inc,_Omega,_omega,_f))
            p = rebound.init_planet(self.sun,a=_a, e=_e, omega=_omega, inc=_inc, 
                                    Omega=_Omega, mass=0., M=_M)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(_a,o.a,places=12)
            self.assertAlmostEqual(_e,o.e,places=12)
            self.assertAlmostEqual(_inc,o.inc,places=12)
            self.assertAlmostEqual(_Omega,rebound.mod2pi(o.Omega),places=12)
            self.assertAlmostEqual(_omega,rebound.mod2pi(o.omega),places=12)
            self.assertAlmostEqual(_M,rebound.mod2pi(o.l-o.Omega-o.omega),places=12)

if __name__ == "__main__":
    unittest.main()