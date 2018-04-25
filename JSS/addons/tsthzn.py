import rebound
import numpy as np


sim = rebound.Simulation()
sim.units = ('AU', 'days', 'Msun')

# We can add Jupiter and four of its moons by name, since REBOUND is linked to the HORIZONS database.
labels = ["Neptune barycenter"]
sim.add(labels, date='2018-03-30 12:00')

sim.status()
print(sim.particles['Neptune barycenter'])


