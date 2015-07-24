REBOUND - An open-source multi-purpose N-body code
==================================================

.. image:: http://img.shields.io/badge/rebound-v2.0.5-green.svg?style=flat
.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat :target: https://github.com/hannorein/rebound/blob/master/LICENSE
.. image:: http://img.shields.io/travis/hannorein/rebound/master.svg?style=flat :target: https://travis-ci.org/hannorein/rebound/
.. image:: http://img.shields.io/badge/arXiv-1110.4876-orange.svg?style=flat :target: http://arxiv.org/abs/1110.4876
.. image:: http://img.shields.io/badge/arXiv-1409.4779-orange.svg?style=flat :target: http://arxiv.org/abs/1409.4779
.. image:: http://img.shields.io/badge/arXiv-1506.01084-orange.svg?style=flat :target: http://arxiv.org/abs/1506.01084


NEW VERSION
-----------
Welcome to REBOUND version 2! We made many changes to the code. Most importanly, REBOUND is now thread-safe and does not use global variables anymore. All the variables that were previously global, are now contained in the `reb_simulation` structure. This has many advantages, for example, you can run separate simulations in parallel from within one process. We also made it possible to choose all modules at runtime (compared to the selection in the Makefile that was used before). This is much more in line with standard UNIX coding practice and does not severely impact the performance (it might even help making REBOUND a tiny bit faster). This makes REBOUND a fully functional shared library. We added a prefix to all public functions and struct definitions: `reb`.

There are still some features that haven't been fully ported. Most importantly, the MPI parallelization and the SWEEP collision detection routine. 

The best way to get and idea of the changes we made is to look at some of the example problems. If you have trouble using the new version or find a bug, please submit an issue or a pull request on github. 

-------------------

.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/shearingsheet.png
.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/dense.png
.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/outersolarsystem.png
.. image:: https://raw.github.com/hannorein/rebound/master/screenshots/disc.png

REBOUND is an N-body integrator, i.e. a software package that can integrate the motion of particles under the influence of gravity. The particles can represent stars, planets, moons, ring or dust particles. REBOUND is very flexible and can be customized to accurately and efficiently solve many problems in astrophysics.  An incomplete feature list of REBOUND:

* Symplectic integrators (WHFast, WH, SEI, LEAPFROG)
* High accuracy non-symplectic integrator with adaptive timestepping (IAS15)
* Support for collisional/granular dynamics, various collision detection routines
* The code is written entirely in C, conforms to the ISO standard C99 and can be used as a thread-safe shared library
* Easy-to-use Python module, installation in 3 words: `pip install rebound`
* Extensive set of example problems in both C and Python
* Real-time, 3D OpenGL visualization (C version)
* Parallelized with OpenMP (for shared memory systems)
* Parallelized with MPI using an essential tree for gravity and collisions (for distributed memory systems)
* No libraries are needed, use of OpenGL/GLUT for visualization is optional
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/rebound
* No configuration is needed to run any of the example problems. Just type `make && ./rebound` in the problem directory to run them
* Comes with standard ASCII or binary output routines 
* Different modules are easily interchangeable at runtime




How to use REBOUND - a quick introduction
-----------------------------------------
    
You can call REBOUND from C or Python. Which programming language you want to use depends on your taste and your specific application. In short: If you simply want to integrate a few particles such as a planetary system with the high order integrator IAS15 or the new symplectic integrator WHFast then use the Python version. If you want to run large simulations with millions of particles, use an exotic integrator, use OpenGL visualizations, or make use of the distributed tree code then use the C version. 

All the computationally expensive parts of REBOUND are written in C. So even if you use the Python version, you'll end up with a very fast code.

To install the *Python version*, simply type the following command into a terminal::

    pip install rebound

To learn more about how to use REBOUND with Python have a look at the iPython/Jupyter tutorials at https://github.com/hannorein/rebound/blob/master/ipython_examples/

To install the *C version*, clone this repository, e.g. by simply copy-and-pasting the following command into your terminal::
    
    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

To learn more about how to use REBOUND with C, study the examples in the `examples/` directory and continue reading this file. You might also want to have a look at the `rebound.h` file in the `src/` directory which contains the API specifications. Last but not least, REBOUND is open source. If you want to know how something works, you can just look at the source code. And of course, you are welcome to e-mail any of the contributors with questions. We'll do our best to answer them quickly.


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

3. Rein & Tamayo (Monthly Notices of the Royal Astronomical Society, Volume 452, Issue 1, p.376-388) describe WHFast, the fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations. http://adsabs.harvard.edu/abs/2015MNRAS.452..376R


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

This section describes the C version of REBOUND. To learn how to install REBOUND for Python have a look at the iPython/Jupiter notebooks at https://github.com/hannorein/rebound/blob/master/ipython_examples/index.ipynb. Hint: It's super easy!

Installation
------------

You can download, compile and run REBOUND on almost any modern operating system within seconds.  Simply copy and paste this line to your terminal and press enter::

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

or if you do not have git installed::

    wget --no-check-certificate https://github.com/hannorein/rebound/tarball/master -O- | tar xvz && cd hannorein-rebound-*/examples/shearing_sheet/ && make && ./rebound

Make sure you have a compiler suite installed. Open a terminal and type `make` and `cc` to test if your installation is complete. If you are on OSX, you can download Xcode from the AppStore (for free). Once installed, open Xcode, go to Settings, then Downloads and install the Command Line Tools. 



Available modules
-----------------

REBOUND is extremely modular. You have the choice between different gravity, collision, boundary and integrator modules. It is also possible to implement completely new modules with minimal effort. In the new version of REBOUND, modules are chosen at runtime by setting flags in the `reb_simulation` structure. 

The following sections list the available modules that come with REBOUND.

**Gravity solvers**::
  
 Module name              | Description
 ------------------------ | -----------
 REB_GRAVITY_COMPENSATED  | Direct summation with compensated summation, O(N^2), default
 REB_GRAVITY_NONE         | No self-gravity
 REB_GRAVITY_BASIC        | Direct summation, O(N^2)
 REB_GRAVITY_TREE         | Oct tree, Barnes & Hut 1986, O(N log(N))
 REB_GRAVITY_OPENCL       | (upgrade to REBOUND 2.0 still in progress) Direct summation, O(N^2), but accelerated using the OpenCL framework.
 REB_GRAVITY_FFT          | (upgrade to REBOUND 2.0 still in progress) Two dimensional gravity solver using FFTW, works in a periodic box and the shearing sheet. 


**Collision detection**::

 Module name              | Description
 ------------------------ | -----------
 REB_COLLISION_NONE       | No collision detection, default
 REB_COLLISION_DIRECT     | Direct nearest neighbour search, O(N^2)
 REB_COLLISION_TREE       | Oct tree, O(N log(N))
 REB_COLLISION_SWEPP      | (upgrade to REBOUND 2.0 still in progress) Plane sweep algorithm, ideal for low dimensional  problems, O(N) or O(N^1.5) depending on geometry 


**Boundaries**::

 Module name              | Description
 ------------------------ | -----------
 REB_BOUNDARY_NONE        | Dummy. Particles are not affected by boundary conditions, default
 REB_BOUNDARY_OPEN        | Particles are removed from the simulation if they leaves the box.
 REB_BOUNDARY_PERIODIC    | Periodic boundary conditions. Particles are reinserted on the other side if they cross the box boundaries. You can use an arbitrary number of ghost-boxes with this module.
 REB_BOUNDARY_SHEAR       | Shear periodic boundary conditions. Similar to periodic boundary conditions, but ghost-boxes are moving with constant speed, set by the shear.
  

**Integrator**::

 Module name              | Description
 ------------------------ | -----------
 REB_INTEGRATOR_IAS15     | IAS15 stands for Integrator with Adaptive Step-size control, 15th order. It is a vey high order, non-symplectic integrator which can handle arbitrary (velocity dependent) forces and is in most cases accurate down to machine precision. IAS15 can integrate variational equations. Rein & Spiegel 2015, Everhart 1985, default
 REB_INTEGRATOR_WHFAST    | WHFast is the integrator described in Rein & Tamayo 2015, it's a second order symplectic Wisdom Holman integrator with 11th order symplectic correctors. It is extremely fast and accurate, uses Gauss f and g functions to solve the Kepler motion and can integrate variational equations.
 REB_INTEGRATOR_EULER     | Euler scheme, first order
 REB_INTEGRATOR_LEAPFROG  | Leap frog, second order, symplectic
 REB_INTEGRATOR_WH        | SWIFT-style Wisdom-Holman Mapping, mixed variable symplectic integrator for the Kepler potential, second order, note that  `integrator_whfast.c` almost always offers better characteristics, Wisdom & Holman 1991, Kinoshita et al 1991
 REB_INTEGRATOR_SEI       | Symplectic Epicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011
 REB_INTEGRATOR_HYBRID    | An experimental hybrid symplectic integrator that uses WHFast for long term integrations but switches over to IAS15 for close encounters.


Code structure
--------------

REBOUND can be used as a shared library. This is *UNIX-way* of using REBOUND. To compile the `librebound.so` file, simply execute `make` in the main directory. However, installing a shared library can sometimes be an obstacle for new users, especially if you want to change the code frequently. For that reason, all the examples that come with REBOUND don't make use of the shared library, but simply compile all the code (including your setup routines) into one single binary file. Here's an example of how to setup a REBOUND simulation:

.. code-block:: c
 
   #include "rebound.h"
   
   int main(int argc, char* argv[]) {
           struct reb_simulation* r = reb_create_simulation();
           r->dt = 0.1;
           r->integrator = REB_INTEGRATOR_WHFAST;
    
           struct reb_particle p1;
           p1.x = 0;  p1.y = 0;  p1.z = 0; 
           p1.vx = 0; p1.vy = 0; p1.vz = 0; 
           p1.m = 1.;
           reb_add(r, p1);
           
           struct reb_particle p2;
           p2.x = 1;  p2.y = 0;  p2.z = 0; 
           p2.vx = 0; p2.vy = 1; p2.vz = 0; 
           p2.m = 0.;
           reb_add(r, p2);
    
           reb_move_to_com(r);    
           reb_integrate(r,100.);
   }

In the first line we include the REBOUND header file. This file contains all the declarationf of the structures and function that we will be using.

Next, we declare the only function in our file. It is the standard C `main()` function. Within that, we first create a `reb_simulation` structure. This is the main structure that contains all the variables, pointers and particles of a REBOUND simulation. You can create multiple `reb_simulation` structures at the same time. The code is thread-safe.

We can then set flags and variables in the `reb_simulation` structure. Note that the `r` variable is a pointer to the structure, so we use the arrow syntax `r->dt = 0.1` to set the variable. The next line chooses the integrator module. Here, we use the WHFast symplectic integrator.
 
We then create two particles, which are represented by the `reb_particle` structure. We set the initial conditions and then add the particle to the simulation using the `reb_add()` function. Note that this function takes two arguments, the first one is the simulation to which you want to add the particle, and the second is the particle that you want to add. 

Finally, we call the REBOUND function `reb_move_to_com()` which moved the particles to a centre of mass reference frame (this prevents particles from drifting away from the origin) and then start the integration. Here, we integrate for 100 time units.

Note that all REBOUND functions start with the three character prefix `reb`. 

Next, let's add a call-back function to the above example. This function will be called after every timestep and we can use it to output simulation data. The relevant function pointer is called `heartbeat` in the `reb_simulation` structure. We first declare and implement the function and then set the pointer in the main routine:

.. code-block:: c

    void heartbeat(struct reb_simulation* r){
           printf("%f\n",r->t);
    }
    int main(int argc, char* argv[]) {
           ...
           r->heartbeat = heartbeat;
           ...
    }

As you can probably guess, this will make the program print out the current time after every timestep. Since the heartbeat function receives the `reb_simulation` structure, you have access to all the variables and particles within the simulation. You don't need any global variables for that. For example, if we wanted to print out the `x` coordinate of the 2nd particle (the index starts at 0, so the second particle has index 1), we could use this heartbeat function.

.. code-block:: c

    void heartbeat(struct reb_simulation* r){
           double x = r->particles[1].x;
           printf("%f\n",x);
    }

REBOUND comes with various built-in output functions that make your life easier. It can for example calculate the orbital elements for you or output to a binary file to save space. The examples are the best way to get to know these functions. You can also look at the `rebound.h` file in the `src/` directory to get an glimpse of the available functions.



Compiling and directory structure
---------------------------------

If you look at the examples in the `examples/` directory, you see one `.c` file and one `Makefile`. All the REBOUND code itself is in the `src/` directory. This setup keeps the directory in which you're working in nice and clean. To compile one of the examples, go to the directory and type `make`. Then the following events happen

* The `Makefile` sets up various environment variables. These determine settings like the compiler optimization flags and which libraries are included (see below). 
* It then creates a symbolic link in the `src/` directory to the `.c` file in the current directory you're in. 
* Next, it calls the `Makefile` in the `src/` directory and compiles the entire code. Note that it compiles everything everytime you execute the script. 
* Finally it copies the binary executable file into the current directory. It's named `rebound`.

You can execute that file with `./rebound`.
If something goes wrong during the compilation of the examples, it is most likely the visualization module. You can turn it off by deleting the line which contains `OPENGL` in the `Makefile`. Of course, you will not see the visualization in real time anymore. See below on how to install GLUT and fix this issue.

If you want to start working on your own problem, simply copy one of the example directories. Then modify the `.c` file and the `Makefile` according to your specific application.  

The other directories are of interest only if you want to use the Python version of REBOUND. More specifically:

* The `rebound/` directory contains python module source files.
* The `python_examples/` directory contains python example problems.
* The `ipython_examples/` directory contains ipython notebooks with examples and tutorials.


Environment variables
---------------------

The makefile in each problem directory sets various environment variables. These determine the compiler optimization flags, the libraries included and basic code settings.

- `export PROFILING=1`. This enables profiling. You can see how much time is spend in the collision, gravity, integrator and visualization modules. This is useful to get an idea about the computational bottleneck.
- `export QUADRUPOLE=0`. This disables the calculation of quadrupole moments for each cell in the tree. The simulation is faster, but less accurate.
- `export OPENGL=1`. This enables real-time OpenGL visualizations and requires both OpenGL and GLUT libraries to be installed. This should work without any further adjustments on any Mac which has Xcode installed. On Linux both libraries must be installed in `/usr/local/`. You can change the default search paths for libraries in the file `src/Makefile`. 
- `export MPI=0`. This disables parallelization with MPI.
- `export OPENMP=1`. This enables parallelization with OpenMP. The number of threads can be set with an environment variable at runtime, e.g.: `export OMP_NUM_THREADS=8`.
- `export CC=gcc`. This flag can be used to override the default compiler. The default compilers are `gcc` for the sequential and `mpicc` for the parallel version. 
- `export LIB=`. Additional search paths for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 
- `export OPT=-O3`. This sets the additional compiler flag `-O3` and optimizes the code for speed. Additional search paths to header files for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 

When you type make in your problem directory, all of these variables are read and passed on to the makefile in the `src/` directory. The `OPENGL` variable, for example, is used to determine if the OpenGL and GLUT libraries should be included. If the variable is `1` the makefile also sets a pre-compiler macro with `-DOPENGL`. Note that because OPENGL is incompatible with MPI, when MPI is turned on (set to 1), OPENGL is automatically turned off (set to 0) in the main makefile. You rarely should have to work directly with the makefile in the `src/` directory yourself.


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
 l       | Toggle limit to screen refresh rate (50Hz/infinity).  

