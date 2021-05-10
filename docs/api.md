# REBOUND API
These pages describe the REBOUND API.
There are two structures (called *objects* in Python) which you will encounter frequently when working with REBOUND.
The first is the [Simulation structure](simulation.md) which contains all the configuration, status and particle data of one REBOUND simulation.
The second is the [Particle structure](particles.md) which represents one particle in a simulation.

REBOUND is a modular code.
You can combine different [gravity solvers](gravity.md), [collision detection algorithms](collisions.md), [boundary conditions](boundaryconditions.md), and [integration methods](integrators.md). 
Not all combinations make physically sense, and not all combinations are supported.
We describe the different modules and their configuration in this section.



!!! Info
    Because the C and Python versions of REBOUND are very similar, we describe both languages in one documentation. 
    The syntax and examples are provided in both C and Python. 
    Use the tabs to switch between them.
