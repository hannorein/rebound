Changelog
=========

This changelog only includes the most important changes in recent updates. For a full log of all changes, please refer to git.

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
>>>>>>> master


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

