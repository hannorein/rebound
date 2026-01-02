import rebound
import unittest
import numpy as np

# Secular modes for Jupiter. Taken from Laskar (1990).
nu5 = np.array([4.2488163, 28.2206942, 3.0895148, 52.1925732, 27.0613982, 29.3799573, 28.8679427, 27.5734578, 5.4070444, 0.6671228]) # frequency, "/yr
A5 = np.array([44119.0e-6, 15750.0e-6, 1800.0e-6, 516.0e-6, 183.0e-6, 178.0e-6, 107.0e-6, 95.0e-6, 62.0e-6, 58.0e-6]) #amplitude
phi5 = np.array([30.676, 308.112, 121.362, 45.551, 218.696, 217.460, 32.614, 43.733, 116.984, 74.116]) # phase, deg
datasep = 120000.0/365.25*2.0*np.pi # 120000 days in units of year/2pi

class TestFrequencyAnalysis(unittest.TestCase):
    def LA1990(self, type):
        Nsamples = 32768
        nfreq = len(nu5)
        inp = np.zeros(Nsamples*2)
        for i in range(nfreq):
            nu = nu5[i]/1296000.0 # to units of radians/(year/2pi)
            inp[0::2] += A5[i]*np.cos(nu*np.arange(Nsamples)*datasep+phi5[i]/180.0*np.pi)
            inp[1::2] += A5[i]*np.sin(nu*np.arange(Nsamples)*datasep+phi5[i]/180.0*np.pi)


        minfreq = 60.0/1296000.0*datasep
        return rebound.frequency_analysis(inp, type=type, minfreq=-minfreq, maxfreq=minfreq)

    def test_LA1990_type_0(self):
        nu, A, phi = self.LA1990(0)
        nfreq = len(nu)
        for i in range(nfreq):
            nu_error = np.abs(nu[i]*1296000.0/datasep-nu5[i])
            self.assertLess(nu_error, 3e-4)
            A_error = np.abs((A[i]-A5[i])/A5[i])
            self.assertLess(nu_error, 2e-3)
            phi_error = np.abs(phi[i]/np.pi*180.0-phi5[i])
            self.assertLess(phi_error, 5e-1)
    def test_LA1990_type_1(self):
        nu, A, phi = self.LA1990(1)
        nfreq = len(nu)
        for i in range(nfreq):
            nu_error = np.abs(nu[i]*1296000.0/datasep-nu5[i])
            self.assertLess(nu_error, 4e-6)
            A_error = np.abs((A[i]-A5[i])/A5[i])
            self.assertLess(nu_error, 1e-5)
            phi_error = np.abs(phi[i]/np.pi*180.0-phi5[i])
            self.assertLess(phi_error, 6e-3)
    def test_LA1990_type_2(self):
        nu, A, phi = self.LA1990(2)
        nfreq = len(nu)
        for i in range(nfreq):
            nu_error = np.abs(nu[i]*1296000.0/datasep-nu5[i])
            self.assertLess(nu_error, 2e-8)
            A_error = np.abs((A[i]-A5[i])/A5[i])
            self.assertLess(nu_error, 3e-7)
            phi_error = np.abs(phi[i]/np.pi*180.0-phi5[i])
            self.assertLess(phi_error, 5e-5)

if __name__ == "__main__":
    unittest.main()
