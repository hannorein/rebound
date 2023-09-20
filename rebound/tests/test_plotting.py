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
            op = rebound.OrbitPlot(self.sim,periastron=True)
            self.assertIsInstance(op.fig,matplotlib.figure.Figure)
            op = rebound.OrbitPlot(self.sim,periastron=True,color=True,unitlabel="AU")
            self.assertIsInstance(op.fig,matplotlib.figure.Figure)
    
    def test_orbitplot_slices(self):
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            import matplotlib; matplotlib.use("pdf")
            op  = rebound.OrbitPlotSet(self.sim,periastron=True)
            self.assertIsInstance(op.fig,matplotlib.figure.Figure)
            op = rebound.OrbitPlot(self.sim,periastron=True,color=True,orbit_style="solid",unitlabel="AU",xlim=[-1.,1])
            self.assertIsInstance(op.fig,matplotlib.figure.Figure)
   

if __name__ == "__main__":
    unittest.main()
