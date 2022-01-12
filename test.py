import rebound
import numpy as np

sim = rebound.Simulation()
sim.add(m=1)
sim.add(a=1,m=1e-3,e=0.9)
sim.add(a=2.341,m=1e-3,e=0.9)
sim.integrator = "IAS15"

ho = sim.create_ode(length=2)

def derivatives_ho(ode, yDot, y, t):
    yDot[0] = y[1]
    yDot[1] = -0.1*y[0]
def energy_ho(ode):
    return 0.5*ode.y[0]**2 + 0.5*ode.y[1]**2

ho.derivatives = derivatives_ho
ho.y[0] = 1
ho.y[1] = 0

E0 = sim.calculate_energy()
E0_ho = energy_ho(ho)
print(sim._odes_N)

times = np.linspace(0.,10.,100)
energies = np.zeros(len(times))
energies_ho = np.zeros(len(times))
x = np.zeros(len(times))
x_ho = np.zeros(len(times))

for i, t in enumerate(times):
    sim.integrate(t)
    x[i] = sim.particles[1].x
    #x_ho[i] = ho.y[0]
    E1 = sim.calculate_energy()
    energies[i] = np.abs((E1-E0)/E0)
    E1_ho = energy_ho(ho)
    energies_ho[i] = np.abs((E1_ho-E0_ho)/E0_ho)
print(sim._odes_N)




