import rebound
import unittest
import os
from rebound.interruptible_pool import InterruptiblePool
        
def runsim(param):
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(a=param)
    sim.integrate(0.1)
    return sim.particles[1].x

class TestInterruptiblePool(unittest.TestCase):
    def test_nopool(self):
        pool = InterruptiblePool(2)
        params = [1.,1.1]
        res = [runsim(params[0]), runsim(params[1])]
        self.assertAlmostEqual(res[0],0.9950041652780258,delta=1e-15)
        self.assertAlmostEqual(res[1],1.095870355119381,delta=1e-15)

    def test_pool(self):
        pool = InterruptiblePool(2)
        params = [1.,1.1]
        res = pool.map(runsim,params)
        self.assertAlmostEqual(res[0],0.9950041652780258,delta=1e-15)
        self.assertAlmostEqual(res[1],1.095870355119381,delta=1e-15)

if __name__ == "__main__":
    unittest.main()
