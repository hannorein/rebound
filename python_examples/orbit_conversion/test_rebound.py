# This script test the orbit conversion routines of REBOUND.
import sys; sys.path.append('../')
import unittest
import rebound
import random
import math

def almost_equal_wrap_2pi(val1,val2, places):
    diff = val2-val1
    diff2 = diff + 2*math.pi if diff < 0 else diff - 2*math.pi
    return True if min(math.fabs(diff), math.fabs(diff2)) < 10**(-places) else False

class cartesian_to_orbital(unittest.TestCase):
    sun = rebound.Particle( m=1.)
    G = 1.
    N_random_tests = 10 # num of random orbits to test in test_rand_r_to_orb_(f or M)
    
    def test_aew2pi(self):
        places=15
        cases = ((0.,10**(-16), True),
                 (0.,10**(-16), True),
                 (0.,10**(-14), False),
                 (0.,10**(-14), False),
                 (0.1,0.2, False))
        for case in cases:
            self.assertIs(almost_equal_wrap_2pi(case[0], 2*math.pi - case[1], places),case[2], '{}'.format(case))
            self.assertIs(almost_equal_wrap_2pi(2*math.pi - case[1], case[0], places),case[2], '{}'.format(case))
            
    def test_keplers_eq(self):
        '''test Kepler's equation'''
        zero_to_pi = [math.pi/10.*q for q in range(11)]
        pi_to_2pi = [math.pi + math.pi/10.*q for q in range(10)]
        e = [0.,0.7,0.9,0.999]
        for ecc in e:
            # E & M match at 0,pi, so if M is in range 0<M<pi (or pi<M<2pi), so should E
            # Also check that Kepler's equation is satisfied to abs precision of 1e-15
            for M in zero_to_pi:
                E = rebound.eccentricAnomaly(ecc,M)
                err = E-ecc*math.sin(E)-M
                self.assertTrue(0.<= E <= math.pi, 
                                "E & M not in same half-plane: E=%.3f, M=%.3f"%(E,M))
                self.assertAlmostEqual(math.fabs(err),0.,places=14,msg=
                "Kepler's Eq not satisfied: e=%.3f, M=%.3f, abs. error=%.3e"%(ecc,M,err))
            for M in pi_to_2pi:
                E = rebound.eccentricAnomaly(ecc,M)
                err = E-ecc*math.sin(E)-M
                self.assertTrue(math.pi<= E <= 2*math.pi,
                                "E & M not in same half-plane: E=%.3f, M=%.3f"%(E,M))       
                self.assertAlmostEqual(math.fabs(err),0.,places=14,msg=
                "Kepler's Eq not satisfied: e=%.3f, M=%.3f, abs. error=%.3e"%(ecc,M,err))

    def test_r_to_orb_defaults(self):
        cases = ({'m':0., 'primary':self.sun, 'a':1},
                 {'m':0., 'primary':self.sun, 'a':1, 'e':0.01},
                 )

        results = ({'a':1.,'anom':0.,'e':0.,'omega':0.,'inc':0.,'Omega':0.},
                   {'a':1.,'anom':0.,'e':0.01,'omega':0.,'inc':0.,'Omega':0.}
                   )
        places=12
        for ctr,kwargs in enumerate(cases):
            p = rebound.kepler_particle(**kwargs)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(results[ctr]['a'],o.a,places=places,msg='{}'.format(kwargs))
            self.assertAlmostEqual(results[ctr]['e'],o.e,places=places,msg='{}'.format(kwargs))
            self.assertAlmostEqual(results[ctr]['inc'],o.inc,places=places,msg='{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(results[ctr]['Omega'],rebound.mod2pi(o.Omega),places), True, '{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(results[ctr]['omega'],rebound.mod2pi(o.omega),places), True, '{}'.format(kwargs))
           
    def test_specified_r_to_orb(self):
        """test conversion from orbital elements to cartesian and back with specified cases"""
        cases = ({'m':0.,'primary':self.sun,'a':1.,'anom':0.,'e':.01,'omega':0.,'inc':0.,'Omega':0.},
                 {'m':0.,'primary':self.sun,'a':1.,'anom':3.,'e':.999,'omega':3.,'inc':2.,'Omega':3.},
                 {'m':0.,'primary':self.sun,'a':1.,'anom':1.728,'e':.851,'omega':1.287,'inc':1.287,'Omega':5.445},
                 #{'m':0.,'primary':self.sun,'a':42.,'anom':0.,'e':1.e-8,'omega':0.,'inc':1.e-8,'Omega':0.}
                 )    
        places=12
        for kwargs in cases:
            p = rebound.kepler_particle(**kwargs)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(kwargs['a'],o.a,places=places, msg='{}'.format(kwargs))
            self.assertAlmostEqual(kwargs['e'],o.e,places=places, msg='{}'.format(kwargs))
            self.assertAlmostEqual(kwargs['inc'],o.inc,places=places, msg='{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['Omega'],rebound.mod2pi(o.Omega),places), True, '{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['omega'],rebound.mod2pi(o.omega),places), True, '{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['anom'],o.f,places), True, '{}'.format(kwargs))
            
    def test_rand_r_to_orb_f(self):    
        places=12
        for q in range(self.N_random_tests):
            kwargs = {'a':random.uniform(1.,2.),'anom':random.uniform(0,2*math.pi),
                      'e':random.uniform(0.,1.),'omega':random.uniform(0,2*math.pi),
                      'inc':random.uniform(0,math.pi),'Omega':random.uniform(0,2*math.pi)}
            
            p = rebound.kepler_particle(m=0.,primary=self.sun,**kwargs)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(kwargs['a'],o.a,places=12,msg='{}'.format(kwargs))
            self.assertAlmostEqual(kwargs['e'],o.e,places=12,msg='{}'.format(kwargs))
            self.assertAlmostEqual(kwargs['inc'],o.inc,places=12,msg='{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['Omega'],rebound.mod2pi(o.Omega),places), True,'{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['omega'],rebound.mod2pi(o.omega),places), True,'{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['anom'],rebound.mod2pi(o.f),places), True,'{}'.format(kwargs))
            
    def test_rand_r_to_orb_f(self):    
        places=12
        for q in range(self.N_random_tests):
            kwargs = {'a':random.uniform(1.,2.),'anom':random.uniform(0,2*math.pi),
                      'e':random.uniform(0.,1.),'omega':random.uniform(0,2*math.pi),
                      'inc':random.uniform(0,math.pi),'Omega':random.uniform(0,2*math.pi),
                      'MEAN':True}
            
            p = rebound.kepler_particle(m=0.,primary=self.sun,**kwargs)
            o = rebound.p2orbit(p,self.sun)
            self.assertAlmostEqual(kwargs['a'],o.a,places=12,msg='{}'.format(kwargs))
            self.assertAlmostEqual(kwargs['e'],o.e,places=12,msg='{}'.format(kwargs))
            self.assertAlmostEqual(kwargs['inc'],o.inc,places=12,msg='{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['Omega'],rebound.mod2pi(o.Omega),places), True,'{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['omega'],rebound.mod2pi(o.omega),places), True,'{}'.format(kwargs))
            self.assertIs(almost_equal_wrap_2pi(kwargs['anom'],rebound.mod2pi(o.l-o.Omega-o.omega),places), True,'{}'.format(kwargs))
        
if __name__ == "__main__":
    unittest.main()
