import numpy as np
import rebound

class State:
    planets = None

    def __init__(self, planets):
        self.planets = planets

    def get_rv(self, times):
        sim = rebound.Simulation()
        sim.add(m=1.)
        for planet in self.planets:
            sim.add(**planet)
        sim.move_to_com()
        
        rv = np.zeros(len(times))
        for i, t in enumerate(times):
            sim.integrate(t)
            rv[i] = sim.particles[0].vx

        return rv

    def get_rv_plotting(self, Npoints=200, tmax=1.5):
        times = np.linspace(0,tmax,Npoints)
        return times, self.get_rv(times)

