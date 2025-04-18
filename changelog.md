# Changelog

This changelog only includes the most important changes in recent updates. For a full log of all changes, please refer to git.

## Version 4.x

### Version 4.4.8
* Added support for symplectic correctors with barycentric coordinates in WHFast.

### Version 4.4.7
* Added option to disable SSL checks for Horizon queries with `rebound.horizons.SSL_CONTEXT = 'unverified'.`
* Added unit tests.
* Added barycentric coordinates for WHFast.
* Bug fix for when MEGNO is used with adaptive timestepping.
* Added more error messages.
* Fixed various issues in documentation.

### Version 4.4.6
* When initializing particles with "uniform" in python, REBOUND now uses its own `reb_random_uniform()` function. This avoids importing the "random" library and makes results reproducible as the random seed of the simulation is used when generating random numbers.
* More cracefull interrupt handling. REBOUND now stop the integration after the next timestep when CTRL-C is pressed the first time. If CTRL-C s pressed twice, then long loops (during gravity, collision calculations) are terminated immediately. Continuing an integration after one CTRL-C press should be easier with this change as the simulation does not get corrupted.
* Fixed typos in documentation.

### Version 4.4.5
* Version updated to test github workflows

### Version 4.4.4
* Fixed several memory leaks and other memory issues. It is unlikely that any of those bugs did affect an simulation.
* When converting units of a particle, the particle radius is now also converted. 
* Added getter/setters for Pal coordinates to the particle structure in python. Syntax is `sim.particles[1].pal_h`, `sim.particles[1].pal_ix`, etc. 

### Version 4.4.3
* REBOUND now raises ImportError if it detects a size mismatch between the C and python Simulation structures.
* Fixes a bug in the WHFast512 synchronization on non-AVX512 systems.
* Fixes a bug in the SimulationArchive in cases where there are multiple snapshots with t=0.
* Updates TRACE switching condition to match Lu et al (2024).
* TRACE binary file size has been reduced.
* Pericenter passage time is now calculated even if particles are not in a Simulation.

### Version 4.4.2
* Fixed bug in TRACE when adding particles.
* Added WHFast fallback for synchronizing WHFast512 simulations with `N_systems` > 1.
* Output version number used to create Simulationarchive if there is a version mismatch.
* Added C example `simulationarchive_fields` which outputs all fields in a simulationarchive for debugging purposes.

### Version 4.4.1
* Fixed bug in TRACE for FULL PERI modes.

### Version 4.4.0
* Added TRACE integrator. See Lu, Hernandez & Rein (2024) for details on this implementation.

### Version 4.3.2
* No longer clipping particles and orbits in visualization.
* Added a scale to visualization. Hide by pressing `t`.
* Option to take a screenshow manually in png (WebGL) or tga (OpenGL) format by pressing `e`.
* Improved `plane` visualization mode. Now supporting hyperbolic orbits. 
* Fixed a memory leak in `reb_simulation_copy`.

### Version 4.3.1
* Added new `plane` visualization mode for orbits. Press `w` to toggle through available orbit visualization modes.
* Added python interface for screenshot API.
* Fixed an issue where no python exception was raised when a particle was added outside a simulation box.
* Renamed `past_N` to breadcrumbs in visualization module.

### Version 4.3.0
* Take screenshots of WebGL based visualizations using the `reb_simulation_output_screenshot()` function. You need to connect one web browser to the simulation in order to take screenshots. 
* Improved synchronization of visualization and simulation on Windows with mutex.
* Fixes an issue that might lead to NaN values when less than the maximum number of planets are used in WHFast512.

### Version 4.2.0
* It is now possible to programmatically change all aspects of a REBOUND visualization. This can be used to set up default viewing options or to render animations. See the C examples in `animation_solar_system` and `animation_saturn_rings`.
* Reworked matrix operations in visualization routines to follow the Model-View-Projection paradigm.
* Fixed an issues where unit tests would fail because a binary file was not deleted.

### Version 4.1.1
* Fixed python wheels for windows.

### Version 4.1.0
* New visualization feature that allows you to show past particle positions and orbits (keyboard commands p, u, and i).
* After pausing a simulation, you can now advance it by a single timestep by pressing the down arrow or 50 timesteps by pressing the page down key.
* Visualization now supports scroll to zoom.
* Fixed memory leaks when using custom ODEs.
* Fixed broken links in documentations.

### Version 4.0.3
* Default IAS15 timestepping criterion is now `adaptive_mode=2`. See Pham, Rein, and Spiegel (2024) for details. To use the old default timestepping criterion, set `adaptive_mode=1`.
* Fixed a race condition that should improve the responsiveness of web based visualizations.
* Removed the glad dependency from emscripten builds which reduces filesize and improves performance.

### Version 4.0.2
* Fixes an issue where the default Makefiles included white spaces after the SERVER and OPENGL variable definitions. This caused the main Makefile to ignore these settings.
* Added `key_callback` function for customizing user interaction in visualizations.
* Added `simulationarchive_viewer` example.
* Included `-sGL_ENABLE_GET_PROC_ADDRESS` flag that is now needed for the latest version of emscripten.

### Version 4.0.1
* Include missing python packages

### Version 4.0.0
* Major API changes and new features! If you have used a previous version of REBOUND, then you will need to update your code. If you have trouble with the migration, open a GitHub issue!
* Many function and variable names have changed. They now follow a coherent naming convention. See the naming convention section in the documentation for more information.
* New visualization module! Previously, using OpenGL visualization required the GLFW library which led to problems on various operating systems. The new visualization module no longer requires ANY dependencies and is compatible with MacOS, Linux, and Windows. It works by running a local web server to which you can point your browser to. In your web browser, an emscripten compiled version of REBOUND handles the WebGL visualization while constantly updating simulation data over HTTP. You can use ssh and port forwarding to visualize simulations on remote servers. Check out the documentation for more details on this new module.
* OpenGL for all the examples has been turned off by default so that new users don't get stuck at this step. To turn on OPENGL simply change the flag in the Makefile.
* Added emscripten support. All C examples (including those using visualizations) are now automatically compiled with emscripten on readthedocs.org so you can run from within the browser. No download or installation required.
* A race condition in OpenGL visualization has been removed. Visualizations run much smoother.
* `reb_random` functions now callable with `r=NULL`. If `r=NULL` then the time and PID is used as a seed.
* Removed support for Simulationarchives with version 2. Added some additional support for reading corrupt/old archives.
* Fixed memory leak in `reb_simulation_copy`.
* Consistent integer sizes for 32/64bit. This includes padding for `reb_particle` which is stored in the Simulationarchive.


## Version 3.x

### Version 3.28.4
* WHFast512 now support the integration of 2 and 4 planet systems in parallel. Providing a speed up of up to 10x.
* The sqrt7 function used by IAS15 now support a wider range of input arguments.

### Version 3.28.3
* Removed distutils requirement in preparation  for python 3.12.
* Removed rebound.InterruptiblePool as it no longer works with recent python version. Updated examples.
* Added Holmberg example.
* Added `adaptive_mode==3` for IAS15 (Aarseth 1985).

### Version 3.28.2
* Implemented own fmemopen implementation on MacOS. This is mainly to appease conda-forge builds.
* Improved sqrt7 algorithm allows larger convergence interval.

### Version 3.28.1
* Improved support for reading old and corrupted Simulationarchives.
* Renamed `ri_ias15.epsilon_global` to `ri_ias15.adaptive_mode`.
* Added new timestep method for IAS15 `ri_ias15.adaptive_mode = 2`. This is experimental for now. Details to be described in Pham, Rein & Spiegel (in prep).
* Added unit tests to check for fused multiply add instruction (these break reproducibility).
* Added phony target in C Makefile to force rebuilding librebound whenever building examples. 

### Version 3.28.0
- Native Windows support. REBOUND can now be built natively on Windows (without WSL) using the Microsoft Visual Studio Compiler. 
- Python Wheels are now provided for Linux, MacOS, and Windows. This should significantly speed up the installation process on a wide variety of systems.

### Version 3.27.0
* In python, Simulation and Particle objects are now picklable. Just like loading Simulations from a binary file, function pointers will need to be re-set manually after unpickling.
* The difference between simulations can now be printed out in a human readable form. Python syntax: `sim.diff(sim2)`. C syntax: `reb_simulation_diff(sim2, sim1, 1)`.
* Reading Simulationarchives with version < 2 is no longer supported.
* The POSIX function fmemopen() is now required to compile REBOUND. This should not affect many users. However, if you are using macOS, the version needs to be >= 10.13 (this version of macOS, High Sierra, was released in 2017). 
* Internal changes on how Simulationarchives are written. 
* Internal variable names that represent the size of allocated buffers now consistently include the name `N_allocated`.
* The TES (Terrestrial Exoplanet Integrator) has been removed. If you wish to use TES, you will need checkout an earlier version.

### Version 3.26.3
* A few more changes to reduce the number of compiler warnings. This should not affect any calculation.

### Version 3.26.2
* Fixed various signed/unsigned int issues. This should reduce the number of compiler warnings but not affect any calculation.

### Version 3.26.1
* Added support for `AVX512` and `FFP_CONTRACT_OFF` environment variables when using pip to install REBOUND.

### Version 3.26.0
* Added WHFast512 integrator (Javaheri, Rein, Tamayo 2023)

### Version 3.25.1
* Bug fixed that prevented the installation via PyPi

### Version 3.25.0
* MPI parts updated and unit tests added
* Fixed machine independence bug in TES.

### Version 3.24.3
* Updated unit tests so they work on 32bit machines

### Version 3.24.2
* Fixed bug in TES ctypes structure

### Version 3.24.1
* Added CORS proxy for Horizons request in pyodide
* Smoother OpenGL animations when using usleep
* TES calculates orbital period automatically

### Version 3.24.0
* Added support for Simulationarchive larger than 4 GB.
* Updated documentation for Lyapunov characteristic number.

### Version 3.23.5
* Added new units shortcuts (year,years,massist)
* Rearranged some loops and switch statements (doesn't affect floating point numbers). 

### Version 3.23.4
* Added pyproject.toml file

### Version 3.23.3
* Changed the way REBOUND reverses the integration direction when the sign of the timestep is inconsistent with respect to the requested final time.
* Fixes a memory leak when a tree code is used
* Fixes an issue where MERCURIUS was not bit-wise reproducible when safe mode was turned off.

### Version 3.23.2
* Minor changes to the python side of Vec3d to make it more compatible with numpy.

### Version 3.23.1
* Minor changes related to the REBOUND Rotations framework.

### Version 3.23.0
* Added the REBOUND Rotations framework.
* Fixes an issue with showing an incorrect periastron location in OrbitPlot for high mass-ratio systems.
* Adds pre and post timestep calls to the ode framework.

### Version 3.22.0
* OrbitPlot is now a class. Checkout the OrbitPlot.ipynb tutorial. This change allows for interactive plots and much faster updates to existing plots. This is great for rendering animations! 

### Version 3.21.0
* Automatic rescaling of first order variational particles has been added. This will allow you to integrate chaotic systems for longer and obtain a more accruate measure of MEGNO and the Lyapunoc exponent. 
* Added `sim.stop()` / `reb_simulation_stop()` to end an integration from within the heartbeat function.

### Version 3.20.1
* Pal coordinates have been added to the `reb_orbit` struct.

### Version 3.20.0
* A new integrator has been added, the Terrestrial Exoplanet Simulation (TES).

### Version 3.19.10
* Fixes another bug int he BS integrator when additional forces are used.

### Version 3.19.9
* Two bugs fixed in the BS integrator. One was related to unitialized memory and the other to issues when the particle number changed.

### Version 3.19.5
* Workaround for urllib support in pyodide added
* Silent warning when InterruptiblePool is not available

### Version 3.19.4
* InterruptiblePool is optional.
* Fixed an issue that occured when switching integrators while using the Simulationarchive.
* Renamed `srand_seed` to make it user accessible. 

### Version 3.19.3
* Added several examples.
* Changed how pypi is rendering the documentation.

### Version 3.19.2
* Fixes a bug relates to test particles of type 0 in MERCURIUS.

### Version 3.19.1
* Some compilers seem to complain that a constant cannot be initialized from a constant. Fixed this so that REBOUND works on colaboratory.

### Version 3.19.0
* Added a Gragg-Bulirsch-Stoer integrator (short BS for Bulirsch-Stoer). This is an adaptive integrator which uses Richardson extrapolation and the modified midpoint method to obtain solutions to ordinary differential equations. The version in REBOUND is based on the method described in Hairer, Norsett, and Wanner 1993 (see section II.9, page 224ff).
* Added the ability to integrate arbitrary ordinary differential equations with REBOUND. The ODEs can be couple to the N-body simulation. This can be used to simulate spin, tides, and other physical effects. The user-defined ODEs are integrated with the new BS integrator.

### Version 3.18.1
* Various improvements and fixes relates to NASA Horizons: small bodies are retrieved correctly, the dates now work with fractional JD values and dates in the format YYYY-MM-DD HH:MM:SS are now supported. 

### Version 3.18.0
* Fixes an issue in the Simulationarchive that prevented REBOUND from seeing more than one snapshot. This only affected simulations with a large number of particles. 

### Version 3.17.5
* REBOUND will now uses the new HTTP API from NASA Horizons. This is significantly faster than the old telnet version. Thanks to Lukas Winkler for implementing this.

### Version 3.17.4
* REBOUND will now attempt to recover binary files and Simulationarchives which have been corrupted. Simulations can be restarted from corrupt files and in most cases the corrupt files will fix themselves.

### Version 3.17.3
* Allow for Horizon queries with future JD dates.

### Version 3.17.2
* Moved some function declarations to rebound.h. This is a temporary fix for REBOUNDx.

### Version 3.17.1
* Fixed an issue where the simulation struct in python did not match the one in C. This might have lead to unexpected behaviour in rare cases.
* Fixed various typos in the documentation
* MERCURIUS switching functions can now be set from Python. Also inluded more built-in switching functions from Hernandez (2019). 

### Version 3.17.0
* Added new 'reb_simulation_add_fmt()' function. This makes adding particles in C as easy as in python.
* Orbits can now also be initialized using the eccentric anomaly.
* Fixed an issue which prevented one loop in the gravity routine form being parallelized with OpenMP.
* Added a warning message when test particles have finite mass.
* More reliable reading of corrupt Simulationarchive files.

### Version 3.16.0
* MERCURIUS: If encounters only involve test-particles (type 0), then the algorithm is now resetting the coordinates of all massive particles after the encounter step. This only changes the outcome at the machine precision, but it makes the trajectories of massive particles independent of the close encounter history. Thanks to Kat Deck for this feature!
* MERCURIUS: The gravity routine is now $O(0.5 \cdot N^2)$ instead of $O(N^2)$ for non-OPENMP runs. This should lead to a noticable improvement in runtime.

### Version 3.15.0
* Orbital parameters of particles can now be changed in-place. For example: 'sim.particles[1].e += 0.1'.
* Implemented more chatty repr functions for most object. Printing REBOUND objects should now give some useful information. 
* Improved support for adding/removing particle in MERCURIUS during collisions.
* REBOUND now outputs an error message when one is trying to remove a particle with a negative index.
* Small updates to the documentation.
* New ipython example added, showing how to use a python collision resolve function.

### Version 3.14.0
* Due to a bug, WHFast was not thread-safe. It is now.
* Random number generator seed is now stored in the Simulationarchive. 
  This allows you to get reproducible random number even after restarting a simulation.
* Random numbers generated with the `reb_rand_*()` functions were not thread-safe.
  They are thread-safe now. Note that this required an API change. All `reb_rand_*()`
  functions now require the simulation structure as an argument. This is because the
  random number generator seed is now stored in the simulation structure. 

### Version 3.13.2
* Correct handling of test particles in reb_transformations.
* Small bug fixes 

### Version 3.13.1
* WHFast: Fixes multiple issues with testparticles in WHFast. 

### Version 3.13.0
* IAS15: Fixes a bug which leads to a biased energy error in long term integrations with fixed timesteps (see Hernandez and Holman 2020). The old version of IAS15 can still be used for the time being by setting ri_ias15.neworder=0.
* IAS15: Does not take variational particles into account when predicting new timesteps. This should be beneficial during close encounters.
* A few improvements have been made to the Simulationarchives code including a more efficient loading procedure for large datasets.

### Version 3.12.3
* Various small bug fixes 
* Added a new function sim.cite() to automatically generate citations depending on the current simulation settings. 

### Version 3.12.2
* Various bug fixes to MERCURIUS
* Performance increase when using the BASIC Gravity Routine with OpenMP
  
### Version 3.12.1
* Bug fixes to LINE and LINETREE algorithms
  
### Version 3.12.0
* Added LINETREE collision search algorithm. 
  This algorithm uses a tree to check if any two particle trajectories overlapped during the last timestep. This
  should be beneficial in large N, low density situation as it allows for much larger timesteps. A modification of the
  collision resolve routine might be necessary to allow for multiple collisions of the same particle during one timestep.
  This depends on the application and the default is to only allow one collision per timestep.

### Version 3.11.1
* Added support for test particles and first-order variational particles to the Embedded Operator Splitting (EOS).
* BASIC Gravity routine changed from O(N^2) to O(0.5 N^2). This should lead to a speed-up in most cases but will break bit-wise reproducibility from earlier versions as the ordering of floating point operations has changed.

### Version 3.11.0
* This version adds the new Embedded Operator Splitting methods from Rein (2019). See the tutorial in the ipython_examples folder for how to use them. 

### Version 3.10.2
* Updates to OrbitPlot. Includes better layout of plot and some syntax changes. See OrbitPlot documentation for the new syntax.

### Version 3.10.1
* Small syntax changes for SABA integrator family. 
* Includes high order integrators by Blanes et al. (2013).

### Version 3.10.0
* Changes for the new version of REBOUNDx. 

### Version 3.9.0
* Added new high order symplectic integrators from Wisdom et al. (1996) and Laskar & Robutel (2001). The implementation of these integrators are discussed in Rein, Tamayo & Brown (2019). 
* Implemented new bit-wise comparison functions for simulations. Python syntax is simply sim1==sim2. 
* Fixed a bug in IAS15 which prevented a restarted simulation to reproduce the original simulation exactly. 

### Version 3.8.3
* Improves and fixes various issues related to variational equations and MEGNO. 

### Version 3.8.2
* Fixes a bug which resulted in duplicate snapshots in Simulationarchives when restarting simulations.

### Version 3.8.1
*   Syntax change on the python side to create a simulation from a binary file or Simulationarchive:
    
    ```python
    rebound.Simulation.from_file("test.bin") becomes rebound.Simulation("test.bin") 
    rebound.Simulation.from_archive("test.bin",5) becomes rebound.Simulation("test.bin",5) 
    ```

### Version 3.8.0
* The hybrid integrator MERCURIUS has been completely rewritten. It can now much more easily be used in simulations where physical collisions occur. There are no more hidden particle arrays in the background, meaning adding and removing particles can occur in the same way as for other integrators. It also works reliably with any additional forces.
* The old hybrid integrator HERMES has been removed. MERCURIUS should always be equal or better in performance and accuracy.

### Version 3.7.1 
* Added getBezierPaths to Simulationarchive to allow for easy plotting of complicated trajectories. To do this, store a lot of snapshots in the Simulationarchive (several per orbit!). 
* Added functionality to add, subtract, multiply and divide simulations. This might be useful when developing new algorithms, but is most likely not useful for most users.

### Version 3.7.0
* Added a deep copy functionality: reb_simulation_copy() in C, and sim.copy() in python. 
* Refactored WHFast to enable calling only certain substeps. 

### Version 3.6.8
* Added the rhill property to reb_orbit in C and the Orbit and Particle classes in Python. This parameter corresponds to the circular Hill radius of the particle: $ a (m/(3M)^{1/3}$.

### Version 3.6.7
* Fixes an issue related to collisions and the Mercurius integrator that prevented the last_collision property to be updated.

### Version 3.6.6
* New: Fancy plotting routine. Usage: rebound.OrbitPlot(sim, fancy=True)

### Version 3.6.5
* One can now add particles from NASA Horizons using Julian Days. For example: sim.add("Earth", date="JD2458327.500000")

### Version 3.6.4
* Fixes a memory leak when using the old Simulationarchive version. Thanks to Ian Rabago for reporting the issue.

### Version 3.6.2
* Fixes a memory leak in the Simulationarchive read function.

### Version 3.6.1
* Removed function calls to open_memstream and fmemopen which might not work on older Mac OSX versions. This only affects the internals and there are no changes to user interface. 
* Minor bug fixes

### Version 3.6.0
* Simulationarchive Version 2. With the new version of the Simulationarchive file format, you can now create snapshots of your simulations without any restrictions. You can change the number of particles, the timestep, even the integrator used during the integration. REBOUND automatically detects what has changed and only stores the differences in incremental snaphots. This reduces the filesize while keeping the format as flexible as possible. The old Simulationarchive Version 1 is still supported for now but might become deprecated in the future. All examples have been updated. As usual these are as usual good starting points for understanding the functionality and the syntax. 

### Version 3.5.12
* Added REB_COLLISION_LINE. This is a collision detection routine which serves for collisions during the last timestep, assuming that all particles travel along straight lines. This can be useful in cases where not every collision needs to be detected exactly, but the overall collision rate should be reproduced. The algorithm is O(N**2).
* Bug related to N_active and variational particles has been fixed.
* A bug where WHFast might not converge in rare cases involving negative timesteps has been fixed.

### Version 3.5.11
* Changed default collision behaviour from hardsphere bouncing to halting the simulation. An exception is raised when using the python version. In C, you need to check the status flag after integrating the simulation.

### Version 3.5.10
* Refactored OrbitPlot.

### Version 3.5.9
* SIGINT handler added. Allows for garceful exit and keyboard interrupts (even from python).

### Version 3.5.8
* WebGL widget text overlay added.

### Version 3.5.7
* Bug fixes related to WebGL widget and ipywidgets version 6

### Version 3.5.6
* Updated WebGL widget to work with ipywidgets version 7

### Version 3.5.5
* Various fixed for Mercurius

### Version 3.5.4
* Bug fix for N_active=-1 (default)

### Version 3.5.3
* Allow for better parallelization of WHFast with OpenMP.
* Addded example of the Solar System with Testparticles.
* Made simulationarchive_append a public function (might be useful for some hacking projects).

### Version 3.5.2
* Fixes an issue with the WebGL widget.
* Fixes an issue with external forces and MERCURIUS.

### Version 3.5.1
* MERCURIUS is not compatible with binary files and the Simulationarchive.

### Version 3.5.0
* The WHFast integrator now supports Jacobi coordinates (default), democratic heliocentric coordinates and WHDS coordinates. The previously separate WHFastHelio integrator has been removed. The coordinate system can now be changed by simply setting the coordinates flag in the ri_whfast struct.
* Included an experimental new integrator MERCURIUS. This is similar to the hybrid integrator in Mercury but uses WHFast and IAS15. Not ready for production yet.

### Version 3.4.0
* Added a screenshot functionality for the WebGL ipython widget. This lets you take screenshots programmatically which is useful to create movies of simulations. 

### Version 3.3.1
* Removed the march=native compiler flag as it seems to be problematic for some OSX/Sierra compilers.

### Version 3.3.0
* JANUS integrator added. This is a bit-wise reversible high-order symplectic integrator. At this time, it remains experimental. Details about this integrator will be published in an upcoming paper.

### Version 3.2.4
* Changes to the WHFastHelio integrator. This integrator now uses democratic heliocentric coordinates and a Hamiltonian splitted as proposed by Hernandez and Dehnen (2017), WHDS, which splits the Hamiltonian into three parts. It has the advantage that the integrator solves the two body problem exactly. It is not compatible with symplectic correctors, this functionality has been removed for WHFastHelio. For very high accuracy integrations of stable planetary systems, the WHFast integrator in Jacobi coordinated (and potentially symplectic correctors) should be better suited.  

### Version 3.2.3
* Various minor bug fixes. Added pre-timestep modifications for REBOUNDx. 

### Version 3.2.2
* Various minor bug fixes. One related to exact_finish_time=1. 

### Version 3.2.0
* Added real-time interactive 3D visualizations using WebGL for Jupyter notebooks. This is an early release. Not everything might be working yet and new feature will be added to the widget class. To try it out, simply run `sim.widget()` in a Jupyter notebook. Note that you need to have ipywidgets installed and enabled. 
* Minor changes to the Visualization backend. This should not have any consequences for users.


### Version 3.1.1
* Now stores the first characters of the current githash in binary files. This is helpful when trying to restart simulations from a binary file and making sure one uses the same version of REBOUND than in the original run. Currently, the git hash is not automatically compared when reloading a binary file. To view the githash, use e.g. hexdump. The hash appears between the first and second zero character in the first 64 bytes of the file. 

### Version 3.1.0
* Updated visualization. REBOUND now uses a modern version of OpenGL (3.3) that allows for custom shaders and therefore better looking visualizations. However, REBOUND now requires glfw3 to compile the visualization module. If you are on a Mac, then the easiest way to install the glfw3 library is with homebrew: `brew tap homebrew/versions && brew install glfw3`. If you are on Linux, you can install it with your package manager, for example with `sudo apt-get install libglfw3-dev`. 

### Version 3.0.0
* Introducing the Simulationarchive. The Simulationarchive allows for exact (bit-by-bit) reproducibility in N-body simulations and a completely new way of analyzing simulations. See Rein&Tamayo (2017) for details.
* The binary format has changed. Binary files created with an earlier version of REBOUND can not be loaded with this version. However, future binary files will be backwards compatible from this point forward.


## Version 2.x
### Version 2.20.6
* Minor bug fixes in HERMES integrator and some examples.

### Version 2.20.5
* NASA Horizons changed a telnet command. This update implements those changes and restores access to NASA Horizons from within REBOUND.

### Version 2.20.4
* Improvements to the Kepler solver. This is typically only relevant for extremly long simulation (1e11 timesteps or more) and extremely accurate simulation with symplectic correctors and a relative energy error of less than 1e-10.

### Version 2.20.3
* Small changes to HERMES integrator. It now has a Solar Switch Factor SSF to allow for close encounters with the central object. 

### Version 2.20.2
* Added adaptive HSF for HERMES integrator. More documentation and paper to follow. 

### Version 2.20.1
* Added symplectic correctors for WHFastHelio integrator. See Wisdom (2006). 
* Improved accuracy of symplectic corrector coefficients for WHFast and WHFastHelio.

### Version 2.20.0
* Added new WHFastHelio integrator. This integrator uses the WHFast Kepler solver, but uses democratic heliocentric coordinates (WHFast itself uses Jacobi coordinates). Heliocentric coordinates are advantages if planets swap positions. 

### Version 2.19.2
* Changes to how particle hashes are handled.

### Version 2.19.1
* This version removes the old SWIFTER based Wisdom-Holman routine, INTEGRATOR_WH. It wasn't working correctly for a while and the WHFast (INTEGRATOR_WHFAST) should be superior in any possible case we can think of. 

### Version 2.19.0
* Added warning/error message system. This allows warning messages to be shown directly in iPython/python programs, rather than being shown on the console. To hide the warning messages, use a filter, e.g.
.. code::  python
    
   with warnings.catch_warnings(record=True) as w:
       warnings.simplefilter("always")
       # Execute a command which triggers a warning message.
       # The message will not show up.
* Improvements regarding the WHFast logic for hyperbolic orbis. No changes should be noticeable to users.

### Version 2.18.9
* Added the reb_simulation_get_serialized_particle_data function for fast access to particle data via numpy array. The full syntax is explained in the documentation. Here is a short example: 
.. code:: python
   
   import numpy as np
   a = np.zeros((sim.N,3),dtype="float64")
   sim.serialize_particle_data(xyz=a)
   print(a)


### Version 2.18.5
* When loading a simulation from a binary file, REBOUND now checks if the version of the binary file is the same as the current version. 
* When saving a simulation to a binary file, all the auxiliary arrays for IAS15 are now stored. This allows for bit-by-bit reproducibility in simulations that are making use of checkpoints.


### Version 2.18.0
* We replaced the old HYBRID integrator with the new and better HERMES integrator. Details of the HERMES integrator will be explained in an upcoming paper Silburt et al (2016, in prep). 

### Version 2.17.0
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
   


### Version 2.0.0
* We made many changes to the code. Most importantly, REBOUND is now thread-safe and does not use global variables anymore. All the variables that were previously global, are now contained in the ``reb_simulation`` structure. This has many advantages, for example, you can run separate simulations in parallel from within one process.
* We also made it possible to choose all modules at runtime (compared to the selection in the ``Makefile`` that was used before). This is much more in line with standard UNIX coding practice and does not severely impact performance (it might even help making REBOUND a tiny bit faster). This makes REBOUND a fully functional shared library. We added a prefix to all public functions and struct definitions: ``reb_``.
* There are still some features that haven't been fully ported. Most importantly, the MPI parallelization and the SWEEP collision detection routine. 
* The best way to get an idea of the changes we made is to look at some of the example problems and the new REBOUND documentation. If you have trouble using the new version or find a bug, please submit an issue or a pull request on github. 

