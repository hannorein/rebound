# Import the rebound module
import rebound
import os

filename = "simulationarchive.bin"
try:
    sim = rebound.Simulation.from_archive(filename)
    sim.simulationarchive_filename = filename
    print("Restarting from simulation archive. Last snapshot found at t=%.1f"%sim.t)
except:
    # Create new simualtion
    print("Creating new simulation.")
    sim = rebound.Simulation()
    sim.add(m=1) # star
    sim.add(m=1e-3, a=1, e=0.01) # planet 1
    sim.add(m=1e-3, a=2.5, e=0.01) # planet 2
    sim.integrator = "whfast"
    sim.dt = 3.1415*2.*6./365.25 # 6 days in units where G=1
    sim.move_to_com()
    sim.initSimulationArchive(filename, interval=2.*3.1415*1e5)

# Run a very long simulation.
# This can be interrupted at any time and then restarted.
sim.integrate(2.*3.1415*1e10) # 10 Gyr

