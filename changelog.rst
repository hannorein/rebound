Changelog
=========

This changelog only includes the most important changes in recent updates. For a full log of all changes, please refer to git.

Version 3.6.1
--------------
* Removed function calls to open_memstream and fmemopen which might not work on older Mac OSX versions. This only affects the internals and there are no changes to user interface. 
* Minor bug fixes

Version 3.6.0
--------------
* SimulationArchive Version 2. With the new version of the SimulationArchive file format, you can now create snapshots of your simulations without any restrictions. You can change the number of particles, the timestep, even the integrator used during the integration. REBOUND automatically detects what has changed and only stores the differences in incremental snaphots. This reduces the filesize while keeping the format as flexible as possible. The old SimulationArchive Version 1 is still supported for now but might become deprecated in the future. All examples have been updated. As usual these are as usual good starting points for understanding the functionality and the syntax. 

Version 3.5.12
--------------
* Added REB_COLLISION_LINE. This is a collision detection routine which serves for collisions during the last timestep, assuming that all particles travel along straight lines. This can be useful in cases where not every collision needs to be detected exactly, but the overall collision rate should be reproduced. The algorithm is O(N**2).
* Bug related to N_active and variational particles has been fixed.
* A bug where WHFast might not converge in rare cases involving negative timesteps has been fixed.

Version 3.5.11
--------------
* Changed default collision behaviour from hardsphere bouncing to halting the simulation. An exception is raised when using the python version. In C, you need to check the status flag after integrating the simulation.

Version 3.5.10
-------------
* Refactored OrbitPlot.

Version 3.5.9
-------------
* SIGINT handler added. Allows for garceful exit and keyboard interrupts (even from python).

Version 3.5.8
-------------
* WebGL widget text overlay added.

Version 3.5.7
-------------
* Bug fixes related to WebGL widget and ipywidgets version 6

Version 3.5.6
-------------
* Updated WebGL widget to work with ipywidgets version 7

Version 3.5.5
-------------
* Various fixed for Mercurius

Version 3.5.4
-------------
* Bug fix for N_active=-1 (default)

Version 3.5.3
-------------
* Allow for better parallelization of WHFast with OpenMP.
* Addded example of the Solar System with Testparticles.
* Made simulationarchive_append a public function (might be useful for some hacking projects).

Version 3.5.2
-------------
* Fixes an issue with the WebGL widget.
* Fixes an issue with external forces and MERCURIUS.

Version 3.5.1
-------------
* MERCURIUS is not compatible with binary files and the SimulationArchive.

Version 3.5.0
-------------
* The WHFast integrator now supports Jacobi coordinates (default), democratic heliocentric coordinates and WHDS coordinates. The previously separate WHFastHelio integrator has been removed. The coordinate system can now be changed by simply setting the coordinates flag in the ri_whfast struct.
* Included a experimental new integrator MERCURIUS. This is similar to the hybrid integrator in Mercury but uses WHFast and IAS15. Not ready for production yet.

Version 3.4.0
-------------
* Added a screenshot functionality for the WebGL ipython widget. This lets you take screenshots programmatically which is useful to create movies of simulations. 

Version 3.3.1
-------------
* Removed the march=native compiler flag as it seems to be problematic for some OSX/Sierra compilers.

Version 3.3.0
-------------
* JANUS integrator added. This is a bit-wise reversible high-order symplectic integrator. At this time, it remains experimental. Details about this integrator will be published in an upcoming paper.

Version 3.2.4
--------------
* Changes to the WHFastHelio integrator. This integrator now uses democratic heliocentric coordinates and a Hamiltonian splitted as proposed by Hernandez and Dehnen (2017), WHDS, which splits the Hamiltonian into three parts. It has the advantage that the integrator solves the two body problem exactly. It is not compatible with symplectic correctors, this functionality has been removed for WHFastHelio. For very high accuracy integrations of stable planetary systems, the WHFast integrator in Jacobi coordinated (and potentially symplectic correctors) should be better suited.  

Version 3.2.3
--------------
* Various minor bug fixes. Added pre-timestep modifications for REBOUNDx. 

Version 3.2.2
--------------
* Various minor bug fixes. One related to exact_finish_time=1. 

Version 3.2.0
--------------
* Added real-time interactive 3D visualizations using WebGL for Jupyter notebooks. This is an early release. Not everything might be working yet and new feature will be added to the widget class. To try it out, simply run `sim.getWidget()` in a Jupyter notebook. Note that you need to have ipywidgets installed and enabled. 
* Minor changes to the Visualization backend. This should not have any consequences for users.


Version 3.1.1
--------------
* Now stores the first characters of the current githash in binary files. This is helpful when trying to restart simulations from a binary file and making sure one uses the same version of REBOUND than in the original run. Currently, the git hash is not automatically compared when reloading a binary file. To view the githash, use e.g. hexdump. The hash appears between the first and second zero character in the first 64 bytes of the file. 

Version 3.1.0
--------------
* Updated visualization. REBOUND now uses a modern version of OpenGL (3.3) that allows for custom shaders and therefore better looking visualizations. However, REBOUND now requires glfw3 to compile the visualization module. If you are on a Mac, then the easiest way to install the glfw3 library is with homebrew: `brew tap homebrew/versions && brew install glfw3`. If you are on Linux, you can install it with your package manager, for example with `sudo apt-get install libglfw3-dev`. 

Version 3.0.0
--------------
* Introducing the Simulation Archive. The Simulation Archive allows for exact (bit-by-bit) reproducibility in N-body simulations and a completely new way of analyzing simulations. See Rein&Tamayo (2017) for details.
* The binary format has changed. Binary files created with an earlier version of REBOUND can not be loaded with this version. However, future binary files will be backwards compatible from this point forward.


Version 2.20.6
--------------
* Minor bug fixes in HERMES integrator and some examples.

Version 2.20.5
--------------
* NASA Horizons changed a telnet command. This update implements those changes and restores access to NASA Horizons from within REBOUND.

Version 2.20.4
--------------
* Improvements to the Kepler solver. This is typically only relevant for extremly long simulation (1e11 timesteps or more) and extremely accurate simulation with symplectic correctors and a relative energy error of less than 1e-10.

Version 2.20.3
--------------
* Small changes to HERMES integrator. It now has a Solar Switch Factor SSF to allow for close encounters with the central object. 

Version 2.20.2
--------------
* Added adaptive HSF for HERMES integrator. More documentation and paper to follow. 

Version 2.20.1
--------------
* Added symplectic correctors for WHFastHelio integrator. See Wisdom (2006). 
* Improved accuracy of symplectic corrector coefficients for WHFast and WHFastHelio.

Version 2.20.0
--------------
* Added new WHFastHelio integrator. This integrator uses the WHFast Kepler solver, but uses democratic heliocentric coordinates (WHFast itself uses Jacobi coordinates). Heliocentric coordinates are advantages if planets swap positions. 

Version 2.19.2
--------------
* Changes to how particle hashes are handled.

Version 2.19.1
--------------
* This version removes the old SWIFTER based Wisdom-Holman routine, INTEGRATOR_WH. It wasn't working correctly for a while and the WHFast (INTEGRATOR_WHFAST) should be superior in any possible case we can think of. 

Version 2.19.0
--------------
* Added warning/error message system. This allows warning messages to be shown directly in iPython/python programs, rather than being shown on the console. To hide the warning messages, use a filter, e.g.
.. code::  python
    
   with warnings.catch_warnings(record=True) as w:
       warnings.simplefilter("always")
       # Execute a command which triggers a warning message.
       # The message will not show up.
* Improvements regarding the WHFast logic for hyperbolic orbis. No changes should be noticable to users.

Version 2.18.9
--------------
* Added the reb_serialize_particle_data function for fast access to particle data via numpy array. The full syntax is explain in the documentation. Here is a short example: 
.. code:: python
   
   import numpy as np
   a = np.zeros((sim.N,3),dtype="float64")
   sim.serialize_particle_data(xyz=a)
   print(a)


Version 2.18.5
--------------
* When loading a simulation from a binary file, REBOUND now checks if the version of the binary file is the same as the current version. 
* When saving a simulation to a binary file, all the auxiliary arrays for IAS15 are now stored. This allows for bit-by-bit reproducability in simulations that are making use of checkpoints.


Version 2.18.0
--------------
* We replaced the old HYBRID integrator with the new and better HERMES integrator. Details of the HERMES integrator will be explained in an upcoming paper Silburt et al (2016, in prep). 

Version 2.17.0
--------------

* What used to be called ``id`` in the particle structure is now called ``hash``. This can be used to uniquely identify particles in a simulation. In many cases, one can just identify particles by their position in the particle array, e.g. using ``sim.particles[5]``. However, in cases where particles might get reordered in the particle array (e.g. when using a tree code), when particles can merge (by using the ``collision_resolve_merge`` routine), or when particles get added or removed manually.
* The syntax is as follows:
.. code:: python
   
   sim = rebound.Simulation()
   sim.add(m=1)
   sim.add(m=1e-3,a=1)
   # Setting a hash using a string:
   sim.particles[1].hash = "planet1"
   # Finding a particle using a string:
   p = sim.get_particle_by_hash("planet1")
   # Setting a random unique hash:
   sim.particles[1].hash = sim.generate_unique_hash() 
   # Save unique hash to find particle later
   uhash = sim.particles[1].hash
   # Find particle using the hash
   p = sim.get_particle_by_hash(uhash)
   


Version 2.0.0
-------------

* We made many changes to the code. Most importanly, REBOUND is now thread-safe and does not use global variables anymore. All the variables that were previously global, are now contained in the ``reb_simulation`` structure. This has many advantages, for example, you can run separate simulations in parallel from within one process.
* We also made it possible to choose all modules at runtime (compared to the selection in the ``Makefile`` that was used before). This is much more in line with standard UNIX coding practice and does not severely impact performance (it might even help making REBOUND a tiny bit faster). This makes REBOUND a fully functional shared library. We added a prefix to all public functions and struct definitions: ``reb_``.
* There are still some features that haven't been fully ported. Most importantly, the MPI parallelization and the SWEEP collision detection routine. 
* The best way to get and idea of the changes we made is to look at some of the example problems and the new REBOUND documentation. If you have trouble using the new version or find a bug, please submit an issue or a pull request on github. 

