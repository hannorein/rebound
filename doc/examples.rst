Examples
========
We provide a lot of examples for REBOUND. We think examples are the easiest way to learn how to use REBOUND.
This is a list of all examples that come with REBOUND. 
You can find all the source code for these in the `examples` and `ipython_examples` directories. 

The examples are sorted by topic. 
Note that some examples use the C version, other the python version of REBOUND. 
Often, the syntax is very similar and you might want to look at the c examples even if you want to write python code and vice verse.

Planetary systems
-----------------
All examples in this section are related to simulations of planetary systems. 

The following examples are good if you're just starting to use REBOUND.

.. toctree::
   c_example_simplest
   ipython/Churyumov-Gerasimenko
   ipython/WHFast
   c_example_solar_system
   c_example_outer_solar_system
   c_example_kozai
   c_example_eccentric_orbit
   c_example_restricted_threebody
   ipython/Horizons
   ipython/Starman

If you want to capture close encounters and collisions between particles, have a look at the following examples.

.. toctree::
   c_example_closeencounter
   c_example_closeencounter_record
   c_example_closeencounter_hybrid
   ipython/CloseEncounters
   c_example_mergers

The following examples demonstrate how to calculate orbital elements.

.. toctree::
   c_example_orbital_elements
   ipython/OrbitalElements
   
If you are interested in simulating a planetary system in which there are many small particles, have a look at the following examples. 

.. toctree::
   c_example_planetesimal_disk_migration
   ipython/EccentricComets
   ipython/PrimordialEarth

Some more advanced topics are covered in the following examples.

.. toctree::
   ipython/AdvWHFast
   ipython/Resonances_of_Jupiters_moons
   ipython/TransitTimingVariations
   ipython/HyperbolicOrbits

Checkpoints and Simulation Archive
----------------------------------
All examples in this section demonstrate how to store, reload and restart simulations. 
They also show how to use the Simulation Archive. 
It is a framework that allows you to not only restart simulations, but recreate them bit-by-bit in a machine independent way. 
With the simulation archive you can run a simulation first, then think later about how you want to analyze it.

.. toctree::
   ipython/Checkpoints
   c_example_simulationarchive
   c_example_restarting_simulation
   ipython/SimulationArchive
   ipython/SimulationArchiveRestart

Variational Equations and Chaos detectors
-----------------------------------------
REBOUND supports first and second order variational equations. 
Varational equations have several advantages over shadow particles when calculating chaos indicators such as Lyapunov exponents. 
They can also be used in optimization problems. 
For more details on variational equations and the math behind them, have a look at the paper `Rein and Tamayo 2016 <https://arxiv.org/abs/1603.03424>`_ and the following examples.

.. toctree::
   c_example_variational_equations
   ipython/VariationalEquations
   ipython/VariationalEquationsWithChainRule

If you are interested in chaos indicators such as MEGNO, have a look at the following examples.

.. toctree::
   c_example_megno
   ipython/Megno
   ipython/PoincareMap
   ipython/FourierSpectrum


Additional forces
-----------------

REBOUND lets you add additional forces to your simulations.
These can be used to simulate radiation drag, general relativistic corrections, planet migration or any other force that you can come up with. 
The following list has a few simple examples. 
You might also want to check out `REBOUNDx <https://github.com/dtamayo/reboundx>`_ which is an add-on to REBOUND that has many additional forces already implemented, so you don't have to do the work. 

.. toctree::
   c_example_J2
   c_example_dragforce
   c_example_prdrag
   c_example_circumplanetarydust
   c_example_planetary_migration
   ipython/Forces


Granular Dynamics
-----------------
The examples in this section show how to use REBOUND for simulations in which particles are colliding with each other.
Note that this is a different type of simulation than simulations of colliding planets.
Here, particles collide often with each other. In the planet case, they collide with each other very rarely. 

.. toctree::
   c_example_bouncing_balls
   c_example_bouncing_balls_corners
   c_example_bouncing_string
   c_example_spreading_ring
   c_example_granulardynamics


Tree code
---------
REBOUND has a built in Barnes-Hut oct tree for collision detection and gravity. 
This allows you, for example, to simulate large gravitationally interacting systems.

.. toctree::
   c_example_selfgravity_disc
   c_example_selfgravity_disc_mpi


Planetary rings
---------------
The following examples show simulations integrating particles in planetary rings. 
They use values characteristic for Saturn's rings but can also be used for other ring systems.

.. toctree::
   c_example_shearing_sheet
   c_example_shearing_sheet_2
   c_example_shearing_sheet_mpi
   c_example_overstability
   ipython/SaturnsRings
   
Visualization
-------------
REBOUND comes with several built-in visualization tools. Most of the C examples can be run with OPENGL visualization (just edit the make file in the corresponding directory). 
If you use Jupyter notebooks, you can use a webGL widget that interactively runs in your notebook. 
REBOUND also comes with tools that can be used to visualize orbits using the matplotlib library. 
The following examples illustrate the visualization feature for python.

.. toctree::
   ipython/WebGLVisualization
   ipython/OrbitPlot

Removing particles
------------------
The following examples illustrate how to remove particles from simulations, i.e. after a collision.

.. toctree::
   c_example_removing_particles_from_simulation
   ipython/EscapingParticles
   ipython/RemovingParticlesFromSimulation

And finally, the following examples didn't fit anywhere else.
They for example show how to use hashes to identify particles.

.. toctree::
   c_example_selfgravity_plummer
   c_example_uniquely_identifying_particles_with_hashes
   c_example_openmp
   c_example_profiling
   c_example_star_of_david
   ipython/Testparticles
   ipython/UniquelyIdentifyingParticlesWithHashes
   ipython/Units
