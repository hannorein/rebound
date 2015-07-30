REBOUND - An open-source multi-purpose N-body code
==================================================

.. image:: http://img.shields.io/badge/rebound-v2.2.1-green.svg?style=flat
.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat :target: https://github.com/hannorein/rebound/blob/master/LICENSE
.. image:: http://img.shields.io/travis/hannorein/rebound/master.svg?style=flat :target: https://travis-ci.org/hannorein/rebound/
.. image:: http://img.shields.io/badge/arXiv-1110.4876-orange.svg?style=flat :target: http://arxiv.org/abs/1110.4876
.. image:: http://img.shields.io/badge/arXiv-1409.4779-orange.svg?style=flat :target: http://arxiv.org/abs/1409.4779
.. image:: http://img.shields.io/badge/arXiv-1506.01084-orange.svg?style=flat :target: http://arxiv.org/abs/1506.01084


NEW VERSION
-----------
Welcome to REBOUND version 2! We made many changes to the code. Most importanly, REBOUND is now thread-safe and does not use global variables anymore. All the variables that were previously global, are now contained in the `reb_simulation` structure. This has many advantages, for example, you can run separate simulations in parallel from within one process. We also made it possible to choose all modules at runtime (compared to the selection in the `Makefile` that was used before). This is much more in line with standard UNIX coding practice and does not severely impact performance (it might even help making REBOUND a tiny bit faster). This makes REBOUND a fully functional shared library. We added a prefix to all public functions and struct definitions: `reb_`.

There are still some features that haven't been fully ported. Most importantly, the MPI parallelization and the SWEEP collision detection routine. 

The best way to get and idea of the changes we made is to look at some of the example problems. If you have trouble using the new version or find a bug, please submit an issue or a pull request on github. 

We're working on a unified way of documenting the code. Stay tuned (or even better, help us!).

-------------------




Examples
========
The following examples can all be found in the `examples` directory. 
Whatever you plan to do with REBOUND, chances are there is already an example available which you can use as a starting point.

* **Bouncing balls.** 

   This example is a simple test of collision detection methods. 

  Directory: examples/bouncing_balls

* **Bouncing balls at corner.** 

   This example tests collision detection methods across box boundaries. There are four particles, one in each corner. To see the ghost boxes in OpenGL press `g` while the simulation is running. 

  Directory: examples/bouncing_balls_corners

* **A string of solid spheres bouncing** 

   This example tests collision detection methods. The example uses a non-square, rectangular box. 10 particles are placed along a line. All except one of the particles are at rest initially. 

  Directory: examples/bouncing_string

* **Radiation forces on circumplanetary dust** 

   This example shows how to integrate circumplanetary dust particles using the IAS15 integrator. The example sets the function pointer `additional_forces` to a function that describes the radiation forces. The example uses a beta parameter of 0.01. The output is custom too, outputting the semi-major axis of every dust particle relative to the planet. 

  Directory: examples/circumplanetarydust

* **Close Encounter** 

   This example integrates a densely packed planetary system which becomes unstable on a timescale of only a few orbits. The IAS15 integrator with adaptive timestepping is used. This integrator automatically decreases the timestep whenever a close encounter happens. IAS15 is very high order and ideally suited for the detection of these kind of encounters. 

  Directory: examples/closeencounter

* **Close Encounter with hybrid integrator (experimental)** 

   This example integrates a densely packed planetary system which becomes unstable on a timescale of only a few orbits. This is a test case for the HYBRID integrator. 

  Directory: examples/closeencounter_hybrid

* **Detect and record close encounters** 

   This example integrates a densely packed planetary system which becomes unstable on a timescale of only a few orbits. The example is identical to the `close_encounter` sample, except that the collisions are recorded and written to a file. What kind of collisions are recorded can be easily modified. It is also possible to implement some additional physics whenever a collision has been detection (e.g. fragmentation). The collision search is by default a direct search, i.e. O(N^2) but can be changed to a tree by using the `collisions_tree.c` module. 

  Directory: examples/closeencounter_record

* **Velocity dependent drag force** 

   This is a very simple example on how to implement a velocity dependent drag force. The example uses the IAS15 integrator, which is ideally suited to handle non-conservative forces. No gravitational forces or collisions are present. 

  Directory: examples/dragforce

* **Example problem: Kozai.** 

   This example uses the IAS15 integrator to simulate a very eccentric planetary orbit. The integrator automatically adjusts the timestep so that the pericentre passages resolved with high accuracy. 

  Directory: examples/eccentric_orbit

* **Granular dynamics.** 

   This example is about granular dynamics. No gravitational forces are present in this example. Two boundary layers made of particles simulate shearing walls. These walls are heating up the particles, create a dense and cool layer in the middle. 

  Directory: examples/granulardynamics

* **J2 precession** 

   This example presents an implementation of the J2 gravitational moment. The equation of motions are integrated with the 15th order IAS15 integrator. The parameters in this example have been chosen to represent those of Saturn, but one can easily change them or even include higher order terms in the multipole expansion. 

  Directory: examples/J2

* **Kozai cycles** 

   This example uses the IAS15 integrator to simulate a Lidov Kozai cycle of a planet perturbed by a distant star. The integrator automatically adjusts the timestep so that even very high eccentricity encounters are resolved with high accuracy. 

  Directory: examples/kozai

* **The chaos indicator MEGNO.** 

   This example uses the IAS15 or WHFAST integrator to calculate the MEGNO of a two planet system. 

  Directory: examples/megno

* **Colliding and merging planets** 

   This example integrates a densely packed planetary system which becomes unstable on a timescale of only a few orbits. The IAS15 integrator with adaptive timestepping is used. The bodies have a finite size and merge if they collide. Note that the size is unphysically large in this example. 

  Directory: examples/mergers

* **Outer Solar System** 

   This example uses the IAS15 integrator to integrate the outer planets of the solar system. The initial conditions are taken from Applegate et al 1986. Pluto is a test particle. This example is a good starting point for any long term orbit integrations. 

   You probably want to turn off the visualization for any serious runs. Go to the makefile and set `OPENGL=0`. 

   The example also works with the WHFAST symplectic integrator. We turn off safe-mode to allow fast and accurate simulations with the symplectic corrector. If an output is required, you need to call ireb_integrator_synchronize() before accessing the particle structure. 

  Directory: examples/outer_solar_system

* **Overstability in Saturn Rings** 

   A narrow box of Saturn's rings is simulated to study the viscous overstability. Collisions are resolved using the plane-sweep method. 

   It takes about 30 orbits for the overstability to occur. You can speed up the calculation by turning off the visualization. Just press `d` while the simulation is running. Press `d` again to turn it back on. 

   You can change the viewing angle of the camera with your mouse or by pressing the `r` key. 

  Directory: examples/overstability

* **How to use unique ids to identify particles** 

   This example shows how to assign ids to particles, and demonstrates different options for removing particles from the simulation. 

  Directory: examples/particles_ids_and_removal

* **Planetary migration in the GJ876 system** 

   This example applies dissipative forces to two bodies orbiting a central object. The forces are specified in terms of damping timescales for the semi-major axis and eccentricity. This mimics planetary migration in a protostellar disc. The example reproduces the study of Lee & Peale (2002) on the formation of the planetary system GJ876. For a comparison, see figure 4 in their paper. The IAS15 or WHFAST integrators can be used. Note that the forces are velocity dependent. Special thanks goes to Willy Kley for helping me to implement the damping terms as actual forces. 

  Directory: examples/planetary_migration

* **Radiation forces** 

   This example provides an implementation of the Poynting-Robertson effect. The code is using the IAS15 integrator which is ideally suited for this velocity dependent force. 

  Directory: examples/prdrag

* **Profiling the shearing sheet example** 

   This example demonstrates how to use the profiling tool that comes with REBOUND to find out which parts of your code are slow. To turn on this option, simple set `PROFILING=1` in the Makefile. Note that enabeling this option makes REBOUND not thread-safe. 

  Directory: examples/profiling

* **Restarting simulations** 

   This example demonstrates how to restart a simulation using a binary file. A shearing sheet ring simulation is used, but the same method can be applied to any other type of simulation. 

  Directory: examples/restarting_simulation

* **Restricted three body problem.** 

   This example simulates a disk of test particles around a central object, being perturbed by a planet. 

  Directory: examples/restricted_threebody

* **Self-gravitating disc.** 

   A self-gravitating disc is integrated using the leap frog integrator. Collisions are not resolved. 

  Directory: examples/selfgravity_disc

* **A self-gravitating Plummer sphere** 

   A self-gravitating Plummer sphere is integrated using the leap frog integrator. Collisions are not resolved. Note that the fixed timestep might not allow you to resolve individual two-body encounters. An alternative integrator is IAS15 which comes with adaptive timestepping. 

  Directory: examples/selfgravity_plummer

* **Shearing sheet (Hill's approximation)** 

   This example simulates a small patch of Saturn's Rings in shearing sheet coordinates. If you have OpenGL enabled, you'll see one copy of the computational domain. Press `g` to see the ghost boxes which are used to calculate gravity and collisions. Particle properties resemble those found in Saturn's rings. 

  Directory: examples/shearing_sheet

* **Shearing sheet (Akihiko Fujii)** 

   This example is identical to the shearing_sheet example but uses a different algorithm for resolving individual collisions. In some cases, this might give more realistic results. Particle properties resemble those found in Saturn's rings. 

   In this collision resolve method, particles are displaced if they overlap. This example also shows how to implement your own collision routine. This is where one could add fragmentation, or merging of particles. 

  Directory: examples/shearing_sheet_2

* **A very simple test problem** 

   We first create a REBOUND simulation, then we add two particles and integrate the system for 100 time units. 

  Directory: examples/simplest

* **Solar System** 

   This example integrates all planets of the Solar System. The data comes from the NASA HORIZONS system. 

  Directory: examples/solar_system

* **Spreading ring** 

   A narrow ring of collisional particles is spreading. 

  Directory: examples/spreading_ring

* **Star of David** 

   This example uses the IAS15 integrator to integrate the "Star od David", a four body system consisting of two binaries orbiting each other. Note that the time is running backwards, which illustrates that IAS15 can handle both forward and backward in time integrations. The initial conditions are by Robert Vanderbei. 

  Directory: examples/star_of_david


OpenGL keyboard command
-----------------------
You can use the following keyboard commands to alter the OpenGL real-time visualizations.::

 Key     | Function
 --------------------------------------------------
 q       | Quit simulation.
 (space) | Pause simulation.
 d       | Pause real-time visualization (simulation continues).
 s       | Toggle three dimensional spheres (looks better)/points (draws faster)
 g       | Toggle ghost boxes
 r       | Reset view. Press multiple times to change orientation.
 x/X     | Move to a coordinate system centred on a particle (note: does not work if particle array is constantly resorted, i.e. in a tree.)
 c       | Toggle clear screen after each time-step.
 w       | Draw orbits as wires (particle with index 0 is central object).  

