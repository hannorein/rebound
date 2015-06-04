REBOUND - An open-source multi-purpose N-body code
==================================================

.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat :target: https://github.com/hannorein/rebound/blob/master/LICENSE>`_
.. image:: http://img.shields.io/badge/arXiv-1110.4876-orange.svg?style=flat :target: http://arxiv.org/abs/1110.4876
.. image:: http://img.shields.io/badge/arXiv-1409.4779-orange.svg?style=flat :target: http://arxiv.org/abs/1409.4779
.. image:: http://img.shields.io/badge/arXiv-1506.01084-orange.svg?style=flat :target: http://arxiv.org/abs/1506.01084

-------------------

.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/shearingsheet.png
.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/dense.png
.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/outersolarsystem.png
.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/disc.png


How to use REBOUND - a quick introduction
-----------------------------------------
    
You can call REBOUND from C or Python. Which programming language you want to use depends on your taste and your specific application. In short: If you simply want to integrate a few particles such as a planetary system with the high order integrator IAS15 or the new symplectic integrator WHFast then use the Python version. If you want to run large simulations with millions of particles, use an exotic integrator, use OpenGL visualizations, or make use of the distributed tree code then use the C version. 

All the computationally expensive parts of REBOUND are written in C. So even if you use the Python version, you'll end up with a very fast code.

To install the Python version, simply type the following command into a terminal::

    pip install rebound

To learn more about how to use REBOUND with Python have a look at the iPython/Jupyter tutorials at https://github.com/hannorein/rebound/blob/master/python_tutorials/

To install the C version, simply copy and paste the following command into your terminal::
    
    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

To learn more about how to use REBOUND with C, continue reading this file.


Feature list 
------------

An incomplete feature list of REBOUND:

* Several symplectic integrators (WHFast, WH, SEI, LEAPFROG)
* High accuracy non-symplectic integrator with adaptive timestepping (IAS15)
* Support for collisional/granular dynamics, various collision detection routines
* The code is written entirely in C, conforms to the ISO standard C99
* Easy-to-use Python module, installation in 3 words: `pip install rebound`
* Extensive set of example problems in both C and Python
* Real-time, 3D OpenGL visualization (C version)
* Parallelized with OpenMP (for shared memory systems)
* Parallelized with MPI using an essential tree for gravity and collisions (for distributed memory systems)
* No libraries are needed, use of OpenGL/GLUT/libpng for visualization is optional
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/rebound
* No configuration is needed to run any of the example problems. Just type `make && ./rebound` in the problem directory to run them
* Standard ASCII or binary output routines 
* Different modules are easily interchangeable by one line in the Makefile


Contributors
------------
* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>
* Akihiko Fujii, National Astronomical Observatory of Japan/University of Tokyo, Tokyo, <akihiko.fujii@nao.ac.jp>
* Dan Tamayo, University of Toronto, <dtamayo@cita.utoronto.ca>

REBOUND is open source. You are invited to contribute to this project if you are using it. Please contact any of the authors above if you have any questions.


Papers
------

There are three papers describing the functionality of REBOUND. 

1. Rein & Liu (Astronomy and Astrophysics, Volume 537, A128, 2012) describe the code structure and the main feature including the gravity and collision routines for many particle systems. http://adsabs.harvard.edu/abs/2012A%26A...537A.128R 

2. Rein & Spiegel (Monthly Notices of the Royal Astronomical Society, Volume 446, Issue 2, p.1424-1437) describe the versatile high order integrator IAS15 which is now part of REBOUND. http://adsabs.harvard.edu/abs/2015MNRAS.446.1424R

3. Rein & Tamayo (Monthly Notices of the Royal Astronomical Society, in press) describe WHFast, the fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations. http://arxiv.org/abs/1506.01084


License
-------
REBOUND is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

REBOUND is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with REBOUND.  If not, see <http://www.gnu.org/licenses/>.


Acknowledgments
---------------
When you use this code or parts of this code for results presented in a scientific publication, please send us a copy of your paper so that we can keep track of all publications that made use of the code. We would greatly appreciate a citation to Rein and Liu (2012) and an acknowledgment of the form:

"Simulations in this paper made use of the REBOUND code which can be downloaded freely at http://github.com/hannorein/rebound."

If you use the IAS15 integrator, please cite Rein and Spiegel (2015).

If you use the WHFast integrator, please cite Rein and Tamayo (2015).


The C version of REBOUND
========================

This section describes the C version of REBOUND. To learn how to install REBOUND for Python have a look at the iPython/Jupiter notebooks at https://github.com/hannorein/rebound/blob/master/python_tutorials/index.ipynb. Hint: It's super easy!

Installation
------------

You can download, compile and run REBOUND on almost any modern operating system within seconds.  Simply copy and paste this line to your terminal and press enter::

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

or if you do not have git installed::

    wget --no-check-certificate https://github.com/hannorein/rebound/tarball/master -O- | tar xvz && cd hannorein-rebound-*/examples/shearing_sheet/ && make && ./rebound

*Note:* Make sure you have a compiler suite installed. Open a terminal and type `make` and `cc` to test if your installation is complete. If you are on OSX, you can download Xcode from the AppStore (for free). Once installed, open Xcode, go to Settings, then Downloads and install the Command Line Tools. 



Available modules
-----------------

REBOUND is extremely modular. You have the choice between different gravity, collision, boundary modules. It is also possible to implement completely new modules with minimal effort. Modules are chosen by setting up symbolic links in the Makefile. There is no need to run a configure script. For example, the Makefile might create a link `gravity.c` that points to one of the gravity modules, say `gravity_tree.c`. This tells the code to use a tree code to do the gravity calculation.

This setup allows you to work on multiple projects at the same time using different modules. When switching to another problem, nothing has to be set-up and the problem can by compiled by simply typing `make` in the corresponding directory (see below).

The following sections list the available modules that come with REBOUND.

**Gravity**::
  
 Module name        | Description
 ------------------ | -----------
 `gravity_none.c`   | No self-gravity
 `gravity_direct.c` | Direct summation, O(N^2)
 `gravity_opencl.c` | Direct summation, O(N^2), but accelerated using the OpenCL framework.
 `gravity_tree.c`   | Oct tree, Barnes & Hut 1986, O(N log(N))
 `gravity_grape.c`  | GRAPE, hardware accelerated direct summation, Sugimoto et al. 1990
 `gravity_fft.c`    | Two dimensional gravity solver using FFTW, works in a periodic box and the shearing sheet. (Not well tested yet.)


**Collision detection**::

 Module name            | Description
 ---------------------- | -----------
 `collisions_none.c`    | No collision detection
 `collisions_direct.c`  | Direct nearest neighbour search, O(N^2)
 `collisions_tree.c`    | Oct tree, O(N log(N))
 `collisions_sweep.c`   | Plane sweep algorithm, ideal for low dimensional  problems, O(N) or O(N^1.5) depending on geometry 
 `collisions_sweepphi.c`| Plane sweep algorithm along the azimuthal angle, ideal for narrow rings in global simulations, O(N) or O(N 1.5) depending on geometry


**Boundaries**::

 Module name            | Description
 ---------------------- | -----------
 `boundaries_open.c`    | Particles are removed from the simulation if they leaves the box.
 `boundaries_none.c`    | Dummy. Particles are not affected by boundary conditions.
 `boundaries_periodic.c`| Periodic boundary conditions. Particles are reinserted on the other side if they cross the box boundaries. You can use an arbitrary number of ghost-boxes with this module.
 `boundaries_shear.c`   | Shear periodic boundary conditions. Similar to periodic boundary conditions, but ghost-boxes are moving with constant speed, set by the shear.
  

Available integrators
---------------------

The following integrators are available within REBOUND. Since May 2015, the integrator can be changed at runtime. Thus, the integrator appears no longer in the Makefile. To set the integrator, set the `integrator` variable in the `probelm_init()` function (see below) to one of the integrator names (it's a C enum)::

 Integrator name   | Description
 ----------------- | -----------
 IAS15             | IAS15 stands for Integrator with Adaptive Step-size control, 15th order. It is a vey high order, non-symplectic integrator which can handle arbitrary (velocity dependent) forces and is in most cases accurate down to machine precision. IAS15 can integrate variational equations. Rein & Spiegel 2015, Everhart 1985
 WHFAST            | WHFast is the integrator described in Rein & Tamayo 2015, it's a second order symplectic Wisdom Holman integrator with 11th order symplectic correctors. It is extremely fast and accurate, uses Gauss f and g functions to solve the Kepler motion and can integrate variational equations.
 EULER             | Euler scheme, first order
 LEAPFROG          | Leap frog, second order, symplectic
 WH                | SWIFT-style Wisdom-Holman Mapping, mixed variable symplectic integrator for the Kepler potential, second order, note that  `integrator_whfast.c` almost always offers better characteristics, Wisdom & Holman 1991, Kinoshita et al 1991
 SEI               | Symplectic Epicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011
 HYBRID            | An experimental hybrid symplectic integrator that uses WHFast for long term integrations but switches over to IAS15 for close encounters.



Directory structure and compilation
-----------------------------------

In the main directory, you find a sub-directory called `src` which contains the bulk parts of the  source code and a directory called `examples` with various example problems. To compile one of the example, you have to go to that directory, for example:

    cd examples/shearing_sheet/

Then, type

    make

This will do the following things    

* It sets various environment variables. These determine settings like the compiler optimization flags and which libraries are included (see below). 
* It creates symbolic links to the active modules. This allows you to choose from different gravity solvers, boundary conditions and collision solvers. For example, to change the gravity solver from using a tree to direct summation you could change `gravity_tree.c` to `gravity_direct.c`. 
* It creates a symbolic link to the current problem file. Each problem file contains the initial conditions and the output routines for the current problem. You do not need to change any file in `src/` to create a new problem unless you want to do something very special. This keeps the initial conditions and the code itself cleanly separated.
* It compiles the code and copies the binary into the current directory.

If something goes wrong, it is most likely the visualization module. You can turn it off by deleting the line which contains `OPENGL` in the makefile. Of course, you will not see the visualization in real time anymore. See below on how to install GLUT and fix this issue.

If you want to start working on your own problem, simply copy one of the example directories. Then modify `problem.c` and `Makefile` according to your application.  


Running REBOUND
---------------

To run the code, simply type

    ./rebound

A window should open and you will see a simulation running in real time. The problem in the directory `examples/shearing_sheet/` simulates the rings of Saturn and uses a local shearing sheet approximation. Have a look at the other examples as well and you will quickly get an idea of what REBOUND can do. 



Environment variables
---------------------

The makefile in each problem directory sets various environment variables. These determine the compiler optimization flags, the libraries included and basic code settings. Let us look at one of the examples `shearing_sheet` in more detail. 

- `export PROFILING=1`. This enables profiling. You can see how much time is spend in the collision, gravity, integrator and visualization modules. This is useful to get an idea about the computational bottleneck.
- `export QUADRUPOLE=0`. This disables the calculation of quadrupole moments for each cell in the tree. The simulation is faster, but less accurate.
- `export OPENGL=1`. This enables real-time OpenGL visualizations and requires both OpenGL and GLUT libraries to be installed. This should work without any further adjustments on any Mac which has Xcode installed. On Linux both libraries must be installed in `/usr/local/`. You can change the default search paths for libraries in the file `src/Makefile`. 
- `export MPI=0`. This disables parallelization with MPI.
- `export OPENMP=1`. This enables parallelization with OpenMP. The number of threads can be set with an environment variable at runtime, e.g.: `export OMP_NUM_THREADS=8`.
- `export CC=gcc`. This flag can be used to override the default compiler. The default compilers are `gcc` for the sequential and `mpicc` for the parallel version. 
- `export LIB=`. Additional search paths for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 
- `export OPT=-O3`. This sets the additional compiler flag `-O3` and optimizes the code for speed. Additional search paths to header files for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 

When you type make in your problem directory, all of these variables are read and passed on to the makefile in the `src/` directory. The `OPENGL` variable, for example, is used to determine if the OpenGL and GLUT libraries should be included. If the variable is `1` the makefile also sets a pre-compiler macro with `-DOPENGL`. Note that because OPENGL is incompatible with MPI, when MPI is turned on (set to 1), OPENGL is automatically turned off (set to 0) in the main makefile. You rarely should have to work directly with the makefile in the `src/` directory yourself.



User-defined functions in the problem.c file
--------------------------------------------

The problem.c file must contain at least three functions. You do need to implement all of them, but a dummy (doing nothing) is sufficient to successfully link the object files. The following documentation describes what these functions do.


- `void problem_init(int argc, char* argv[])`

    This routine is where you read command line arguments and set up your initial conditions. REBOUND does not come with a built-in functionality to read configuration files at run-time. We consider this not a missing feature. In REBOUND, you have one `problem.c` file for each problem. Thus, everything can be set within this file. There are, of course, situation in which you want to do something like a parameter space survey. In almost all cases, you vary only a few parameters. You can easily read these parameters from the command line.
 
    Here is an example that reads in a command line argument given to REBOUND in the standard unix format `./rebound --boxsize=200.`. A default value of 100 is used if no parameter is passed to REBOUND.::

        // At the top of the problem.c file add
        #include "input.h"
        // In problem_init() add
        boxsize = input_get_double(argc,argv,"boxsize",100.);

- `void problem_output()`

    This function is called at the beginning of the simulation and at the end of each time-step. You can implement your output routines here. Many basic output functions are already implemented in REBOUND. See `output.h` for more details. The function `output_check(odt)` can be used to easily check if an output is needed if you want to trigger and output once per time interval `odt`. For example, the following code snippet outputs some timing statistics to the console every 10 time-steps::
    
        if (output_check(10.*dt)){
            output_timing();
        }
 
- `void problem_finish()`

    This function is called at the end of the simulation, when t >= tmax. This is the last chance to output any quantities before the program ends.


- `void problem_additional_forces()` (optional function pointer)

    In addition to the four mandatory functions that need to be present, you can also define some other functions and make them callable by setting a function pointer. The function pointer `problem_additional_forces()` which is called one or more times per time-step whenever the forces are updated. This is where you can implement all kind of things such as additional forces onto particles. 
    
    The following lines of code implement a simple velocity dependent force.  `IAS15` is best suited for this (see `examples/dragforce`)::
    
        void velocity_dependent_force(){
            for (int i=1;i<N;i++){
               particles[i].ax -= 0.0000001 * particles[i].vx;
               particles[i].ay -= 0.0000001 * particles[i].vy;
               particles[i].az -= 0.0000001 * particles[i].vz;
            }
        }
    
    Make sure you set the function pointer in the `problem_init()` routine::
    
        problem_additional_forces = velocity_dependent_force;
    
    By default, all integrators assume that the forces are velocity dependent. If all forces acting on particles only depend on positions, you can set the following variable (defined in `integrator.h`) to `0` to speed up the calculation::
    
        // Add to problem_init()
        integrator_force_is_velocitydependent = 0;


How to install GLUT 
-------------------

The OpenGL Utility Toolkit (GLUT) comes pre-installed as a framework on Mac OSX. If you are working on another operating system, you might have to install GLUT yourself if you see an error message such as `error: GL/glut.h: No such file or directory`. On Debian and Ubuntu, simply make sure the `freeglut3-dev` package is installed. If glut is not available in your package manager, go to http://freeglut.sourceforge.net/ download the latest version, configure it with `./configure` and compile it with `make`. Finally install the library and header files with `make install`. 

You can also install freeglut in a non-default installation directory if you do not have super-user rights by running the freeglut installation script with the prefix option::

    mkdir ${HOME}/local
    ./configure --prefix=${HOME}/local
    make all && make install

Then, add the following lines to the REBOUND Makefile::

    OPT += -I$(HOME)/local/include
    LIB += -L$(HOME)/local/lib

Note that you can still compile and run REBOUND even if you do not have GLUT installed. Simply set `OPENGL=0` in the makefile (see below). 


Examples
========
The following examples can all be found in the `examples` directory. 
Whatever you plan to do with REBOUND, chances are there is already an example available which you can use as a starting point.


examples/bouncing_balls
  This example is a simple test of collision detection
  methods. To change the collision detection algorithm, you can replace
  the module collisions_direct.c to either collisions_tree.c or
  collisions_sweep.c in the Makefile.
  
  Modules used: ``gravity_direct.c`` ``boundaries_periodic.c`` ``collisions_direct.c``.

examples/bouncing_balls_corners
  This example tests collision detection methods across box boundaries.
  There are four particles, one in each corner. To see the ghost boxes in OpenGL
  press `g` while the simulation is running.
  
  Modules used: ``gravity_direct.c`` ``boundaries_periodic.c`` ``collisions_tree.c``.

examples/bouncing_string
  This example tests collision detection methods.
  The example uses a non-square, rectangular box. 10 particles are placed
  along a line. All except one of the particles are at rest
  initially.
  
  Modules used: ``gravity_none.c`` ``boundaries_periodic.c`` ``collisions_direct.c``.

examples/circumplanetarydust
  This example shows how to integrate circumplanetary
  dust particles using the `integrator_ias15.c` module.
  The example sets the function pointer `problem_additional_forces`
  to its own function that describes the radiation forces.
  The example uses a beta parameter of 0.01.
  The output is custom too, outputting the semi-major axis of
  every dust particle relative to the planet.
  Only one dust particle is used in this example, but there could be
  many.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/closeencounter
  This example integrates a densely packed planetary system
  which becomes unstable on a timescale of only a few orbits. The IAS15
  integrator with adaptive timestepping is used. This integrator
  automatically decreases the timestep whenever a close
  encounter happens. IAS15 is very high order and ideally suited for the
  detection of these kind of encounters.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/closeencounter_hybrid
  This example integrates a densely packed planetary system
  which becomes unstable on a timescale of only a few orbits.
  This is a test case for the HYBRID integrator.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/closeencounter_record
  This example integrates a densely packed planetary system
  which becomes unstable on a timescale of only a few orbits.
  The example is identical to the `close_encounter` sample, except that
  the collisions are recorded and written to a file. What kind of collisions
  are recorded can be easily modified. It is also possible to implement some
  additional physics whenever a collision has been detection (e.g. fragmentation).
  The collision search is by default a direct search, i.e. O(N^2) but can be
  changed to a tree by using the `collisions_tree.c` module.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_direct.c``.

examples/dragforce
  This is a very simple example on how to implement a velocity
  dependent drag force. The example uses the IAS15 integrator, which
  is ideally suited to handle non-conservative forces.
  No gravitational forces or collisions are present.
  
  Modules used: ``gravity_none.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/eccentric_orbit
  This example uses the IAS15 integrator to simulate
  a very eccentric planetary orbit. The integrator
  automatically adjusts the timestep so that the pericentre passages
  resolved with high accuracy.
  
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/forced_migration
  This example applies dissipative forces to two
  bodies orbiting a central object. The forces are specified
  in terms of damping timescales for the semi-major axis and
  eccentricity. This mimics planetary migration in a protostellar disc.
  The example reproduces the study of Lee & Peale (2002) on the
  formation of the planetary system GJ876. For a comparison,
  see figure 4 in their paper. The IAS15 integrator is used
  because the forces are velocity dependent.
  Special thanks goes to Willy Kley for helping me to implement
  the damping terms as actual forces.
  
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/granulardynamics
  This example is about granular dynamics. No gravitational
  forces are present in this example, which is why the module
  `gravity_none.c` is used. Two boundary layers made of
  particles simulate shearing walls. These walls are heating
  up the particles, create a dense and cool layer in the middle.
  
  Modules used: ``gravity_none.c`` ``boundaries_periodic.c`` ``collisions_tree.c``.

examples/J2
  This example presents an implementation of the J2
  gravitational moment. The equation of motions are integrated with
  the 15th order IAS15 integrator. The parameters in this examples
  have been chosen to represent those of Saturn, but you can easily
  change them or even include higher order terms in the multipole
  expansion.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/kozai
  This example uses the IAS15 integrator to simulate
  a Lidov Kozai cycle of a planet perturbed by a distant star. The integrator
  automatically adjusts the timestep so that even very high
  eccentricity encounters are resolved with high accuracy.
  
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/megno
  This example uses the IAS15 integrator
  to calculate the MEGNO of a two planet system.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/mergers
  This example integrates a densely packed planetary system
  which becomes unstable on a timescale of only a few orbits. The IAS15
  integrator with adaptive timestepping is used. The bodies have a finite
  size and merge if they collide. Note that the size is unphysically large
  in this example.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_direct.c``.

examples/opencl
  A self-gravitating disc is integrated using
  the OpenCL direct gravity summation module.
  
  This is a very simple implementation (see `gravity_opencl.c`).
  Currently it only supports floating point precision. It also
  transfers the data back and forth from the GPU every timestep.
  There are considerable improvements to be made. This is just a
  proof of concept. Also note that the code required N to be a
  multiple of the workgroup size.
  
  You can test the performance increase by running:
  `make direct && ./rebound`, which will run on the CPU and
  `make && ./rebound`, which will run on the GPU.
  
  The Makefile is working with the Apple LLVM compiler. Changes
  might be necessary for other compilers such as gcc.
  
  
  Modules used: ``gravity_opencl.c`` ``boundaries_open.c`` ``collisions_none.c`` ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/outer_solar_system
  This example uses the IAS15 integrator
  to integrate the outer planets of the solar system. The initial
  conditions are taken from Applegate et al 1986. Pluto is a test
  particle. This example is a good starting point for any long term orbit
  integrations.
  
  You probably want to turn off the visualization for any serious runs.
  Just go to the makefile and set `OPENGL=0`.
  
  The example also works with the Wisdom-Holman symplectic integrator.
  Simply change the integrator to `integrator_wh.c` in the Makefile.
  
  Modules used: ``gravity_direct.c`` ``boundaries_none.c`` ``collisions_none.c``.

examples/overstability
  A narrow box of Saturn's rings is simulated to
  study the viscous overstability. Collisions are resolved using
  the plane-sweep method.
  
  It takes about 30 orbits for the overstability to occur. You can
  speed up the calculation by turning off the visualization. Just press
  `d` while the simulation is running. Press `d` again to turn it back on.
  
  You can change the viewing angle of the camera with your mouse or by pressing
  the `r` key.
  
  Modules used: ``gravity_none.c`` ``boundaries_shear.c`` ``collisions_sweep.c``.

examples/prdrag
  This example provides an implementation of the
  Poynting-Robertson effect. The code is using the IAS15 integrator
  which is ideally suited for this velocity dependent force.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/restarting_simulation
  This example demonstrates how to restart a simulation
  using a binary file. A shearing sheet ring simulation is used, but
  the same method can be applied to any other type of simulation.
  
  First, run the program with `./rebound`.
  Random initial conditions are created and
  a restart file is written once per orbit.
  Then, to restart the simulation, run the
  program with `./rebound --restart restart.bin`.
  
  
  Modules used: ``gravity_direct.c`` ``boundaries_shear.c`` ``collisions_direct.c``.

examples/restricted_threebody
  This example simulates a disk of test particles around
  a central object, being perturbed by a planet.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/restricted_threebody_mpi
  This problem uses MPI to calculate the restricted three
  body problem. Active particles are copied to all nodes. All other
  particles only exist on one node and are not automatically (re-)
  distributed. There is not domain decomposition used in this example.
  Run with `mpirun -np 4 nbody`.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/selfgravity_disc
  A self-gravitating disc is integrated using
  the leap frog integrator. Collisions are not resolved.
  
  Modules used: ``gravity_tree.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/selfgravity_disc_grape
  A self-gravitating disc is integrated using
  the leap frog integrator. This example is using the GRAPE
  module to calculate the self-gravity. You need to have a physical
  GRAPE card in your computer to run this example.
  Collisions are not resolved.
  
  Modules used: ``gravity_grape.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/selfgravity_plummer
  A self-gravitating Plummer sphere is integrated using
  the leap frog integrator. Collisions are not resolved. Note that the
  fixed timestep might not allow you to resolve individual two-body
  encounters. An alternative integrator is `integrator_ias15.c` which
  comes with adaptive timestepping.
  
  Modules used: ``gravity_tree.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/shearing_sheet
  This example simulates a small patch of Saturn's
  Rings in shearing sheet coordinates. If you have OpenGL enabled,
  you'll see one copy of the computational domain. Press `g` to see
  the ghost boxes which are used to calculate gravity and collisions.
  Particle properties resemble those found in Saturn's rings.
  
  
  Modules used: ``gravity_tree.c`` ``boundaries_shear.c`` ``collisions_tree.c``.

examples/shearing_sheet_2
  This example is identical to the shearing_sheet
  example but uses a different algorithm for resolving individual
  collisions. In some cases, this might give more realistic results.
  Particle properties resemble those found in Saturn's rings.
  
  In this collision resolve method, particles are displaced if they
  overlap. This example also shows how to implement your own collision
  routine. This is where one could add fragmentation, or merging of
  particles.
  
  
  Modules used: ``gravity_tree.c`` ``boundaries_shear.c`` ``collisions_tree.c``.

examples/shearing_sheet_fft
  This problem is identical to the other shearing
  sheet examples but uses an FFT based gravity solver.
  To run this example, you need to install the FFTW library.
  Collisions are detected using a plane sweep algorithm.
  There is no tree present in this simulation.
  
  Modules used: ``gravity_fft.c`` ``boundaries_shear.c`` ``collisions_sweep.c``.

examples/shearing_sheet_grape
  This is yet another shearing sheet example,
  it uses a GRAPE to calculate gravity. Note that you need to have
  a physical GRAPE card installed in your computer to run this
  simulation. Particle properties resemble those found in
  Saturn's rings.
  
  Modules used: ``gravity_grape.c`` ``boundaries_shear.c`` ``collisions_sweep.c``.

examples/shearing_sheet_profiling
  This example demonstrates how to use the
  profiling tool that comes with REBOUND to find out which parts
  of your code are slow. To turn on this option, simple set
  `PROFILING=1` in the Makefile.
  
  Modules used: ``gravity_tree.c`` ``boundaries_shear.c`` ``collisions_tree.c``.

examples/simple
  This example uses the IAS15 integrator
  to calculate the MEGNO of a two planet system.
  
  Modules used: ``gravity_direct.c`` ``boundaries_none.c`` ``collisions_none.c``.

examples/solar_system
  This example integrates all planets of the Solar
  System. The data comes from the NASA HORIZONS system.
  
  Modules used: ``gravity_direct.c`` ``boundaries_none.c`` ``collisions_none.c``.

examples/spreading_ring
  A narrow ring of collisional particles is spreading.
  The example uses the Wisdom Holman integrator. A plane-sweep algorithm
  in the phi direction is used to detect collisions.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_sweepphi.c``.

examples/star_of_david
  This example uses the IAS15 integrator
  to integrate the "Star od David", a four body system consisting of two
  binaries orbiting each other. Note that the time is running backwards,
  which illustrates that IAS15 can handle both forward and backward in time
  integrations. The initial conditions are by Robert Vanderbei. For more
  information see http://www.princeton.edu/%7Ervdb/WebGL/New.html
  
  Modules used: ``gravity_direct.c`` ``boundaries_none.c`` ``collisions_none.c``.

examples/stark
  This example calculates the Stark problem.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/stark_radial
  This example uses the IAS15 integrator
  to calculate the MEGNO of a two planet system.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/symplectic_integrator
  This example uses the symplectic Wisdom Holman (WH) integrator
  to integrate test particles on eccentric orbits in a fixed potential.
  Note that the WH integrator assumes that the central object is at the origin.
  
  Modules used: ``gravity_direct.c`` ``boundaries_open.c`` ``collisions_none.c``.

examples/viewer
  This example doesn't simulate anything. It's just a
  visualization toll that can display data in the form x, y, z, r.
  This might be useful when large simulations have been run and you want
  to look (at parts of) it at a later time.
  
  Note that this example uses only dummy modules.
  
  Modules used: ``gravity_none.c`` ``boundaries_periodic.c`` ``collisions_dummy.c``.

examples/whfast
  This example uses the symplectic Wisdom Holman (WH) integrator
  to integrate test particles on eccentric orbits in a fixed potential.
  Note that the WH integrator assumes that the central object is at the origin.
  
  Modules used: ``gravity_direct.c`` ``boundaries_none.c`` ``collisions_none.c``.

OpenGL keyboard command
-----------------------
You can use the following keyboard commands to alter the OpenGL real-time visualizations.::

 Key     | Function
 -------------------------
 (space) | Pause simulation.
 d       | Pause real-time visualization (simulation continues).
 q       | Quit simulation.
 s       | Toggle three dimensional spheres (looks better)/points (draws faster)
 g       | Toggle ghost boxes
 r       | Reset view. Press multiple times to change orientation.
 x/X     | Move to a coordinate system centred on a particle (note: does not work if particle array is constantly resorted, i.e. in a tree.)
 t       | Show tree structure.
 m       | Show centre of mass in tree structure (only available when t is toggled on).
 p       | Save screen shot to file.
 c       | Toggle clear screen after each time-step.
 w       | Draw orbits as wires (particle with index 0 is central object).  

