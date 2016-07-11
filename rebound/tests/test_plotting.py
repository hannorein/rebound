import rebound
import unittest
import warnings

class TestPlotting(unittest.TestCase):
    def setUp(self):
        self.sim = rebound.Simulation()
        self.sim.add(m=1)
        self.sim.add(m=1e-3,a=1,e=0.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
        self.sim.add(m=1e-3,a=-2,e=1.1,omega=0.1,M=0.1,inc=0.1,Omega=0.1)
    
    def tearDown(self):
        self.sim = None
    
    def test_orbitplot(self):
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            import matplotlib; matplotlib.use("pdf")
            import numpy as np
            t = np.array(1.)
            plot = rebound.OrbitPlot(self.sim,periastron=True)
            self.assertIsInstance(plot,matplotlib.figure.Figure)
            plot = rebound.OrbitPlot(self.sim,periastron=True,color=True,trails=True,unitlabel="AU")
            self.assertIsInstance(plot,matplotlib.figure.Figure)
    
    def test_orbitplot_slices(self):
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            import matplotlib; matplotlib.use("pdf")
            import numpy as np
            t = np.array(1.)
            plot = rebound.OrbitPlot(self.sim,periastron=True,slices=True)
            self.assertIsInstance(plot,matplotlib.figure.Figure)
            plot = rebound.OrbitPlot(self.sim,periastron=True,color=True,trails=True,unitlabel="AU",slices=True,limz=1.)
            self.assertIsInstance(plot,matplotlib.figure.Figure)
   

if __name__ == "__main__":
    unittest.main()
