import numpy as np
import rebound

class State:
    planets = None

    def __init__(self, planets):
        self.planets = planets

    def setup_sim(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        for planet in self.planets:
            sim.add(**planet)
        sim.move_to_com()
        return sim

    def get_rv(self, times):
        sim = self.setup_sim()
        
        rv = np.zeros(len(times))
        for i, t in enumerate(times):
            sim.integrate(t)
            rv[i] = sim.particles[0].vx

        return rv

    def get_rv_plotting(self, Npoints=200, tmax=1.5):
        times = np.linspace(0,tmax,Npoints)
        return times, self.get_rv(times)


    def get_chi2(self, obs):
        sim = self.setup_sim()

        chi2 = 0.
        for i, t in enumerate(obs.t):
            sim.integrate(t)
            chi2 += (sim.particles[0].vx-obs.rv[i])**2
        return chi2/(obs.error**2 * obs.Npoints)



