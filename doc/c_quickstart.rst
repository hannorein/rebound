Quick User Quide (C)
====================

This section describes the C version of REBOUND. 

Installation
------------

You can download, compile and run REBOUND on almost any modern operating system within seconds.  Simply copy and paste this line to your terminal and press enter::

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

or if you do not have git installed::

    wget --no-check-certificate https://github.com/hannorein/rebound/tarball/master -O- | tar xvz && cd hannorein-rebound-*/examples/shearing_sheet/ && make && ./rebound

Make sure you have a compiler suite installed. Open a terminal and type `make` and `cc` to test if your installation is complete. If you are on OSX, you can download Xcode from the AppStore (for free). Once installed, open Xcode, go to Settings, then Downloads and install the Command Line Tools. 

Note:  REBOUND does not work on Windows, and we currently do not have plans to support it.

Code structure
--------------

REBOUND can be used as a shared library. However, installing a system-wide shared library can sometimes be an obstacle for new users, especially if you want to change the code frequently or don't have root access. For that reason, all the examples can be compiled by simply typing `make` in any of the example directories.

Let's look at how to setup a simple REBOUND simulation:

.. code-block:: c
 
   #include "rebound.h"
   
   int main(int argc, char* argv[]) {
           struct reb_simulation* r = reb_create_simulation();
           r->dt = 0.1;
           r->integrator = REB_INTEGRATOR_WHFAST;
    
           struct reb_particle p1 = {0};
           p1.m = 1.;
           reb_add(r, p1);
           
           struct reb_particle p2 = {0};
           p2.x = 1;
           p2.vy = 1;
           p2.m = 0.;
           reb_add(r, p2);
    
           reb_move_to_com(r);    
           reb_integrate(r,100.);
   }

In the first line we include the REBOUND header file. This file contains all the declarations of the structures and functions that we will be using.

Next, we declare the only function in our file. It is the standard C `main()` function. Within that, we first create a `reb_simulation` structure. This is the main structure that contains all the variables, pointers and particles of a REBOUND simulation. You can create multiple `reb_simulation` structures at the same time. REBOUND is thread-safe.

We can then set flags and variables in the `reb_simulation` structure. Note that the `r` variable is a pointer to the structure, so we use the arrow syntax `r->` to set variables. The next line chooses the integrator module. Here, we use the WHFast symplectic integrator.
 
We then create two particles, both of which are represented by a `reb_particle` structure. The `= {0}` syntax ensures that our structs are initialized with zeros. We set the initial conditions (the ones we don't want to be zero) and then add the particle to the simulation using the `reb_add()` function. Note that this function takes two arguments, the first one is the simulation to which you want to add the particle, and the second is the particle that you want to add. 

Finally, we call the REBOUND function `reb_move_to_com()`. It moves the particles to a centre of mass reference frame (this prevents particles from drifting away from the origin). We then start the integration. Here, we integrate for 100 time units. By default REBOUND used units in which G=1, thus a particle around an m=1 mass central object at a semi-major axis of 1 needs 2pi time units for one orbit.

Note that all REBOUND functions start with the three character prefix `reb_`. 

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
* Next, it calls the `Makefile` in the `src/` directory and compiles the entire REBOUND code into a shared library. 
* It then creates a symbolic link from the current directory to the location of the share library in the src directory. 
* Finally it compiles your code, the `problem.c` file, into an executable file. 

You can execute that file with `./rebound`.
After you edited a file, you can simply type `make` again to recompile.
If you change any of the environment variables, clean the build directiory first, by executing `make clean`.

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

To rotate the view, simple trag the simulation with the mouse. To zoom in, press the shift key and then trag the simulation with the mouse. 


API Documentation
-----------------
We provide a full API documentation in a separate file. The most important REBOUND API structures and functions are listed below. 
Note that you can also look at the code itself. The starting point is the `rebound.h` file in the `src/` directory. 
This is where the public API is defined. 

The reb_simulation structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenstruct:: reb_simulation
   :members:

Main REBOUND functions
^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: MainRebFunctions
   :members:

Tool functions
^^^^^^^^^^^^^^

.. doxygengroup:: ToolsRebFunctions
   :members:

Output functions
^^^^^^^^^^^^^^^^

.. doxygengroup:: OutputRebFunctions
   :members:

Particle setup functions
^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: SetupRebFunctions
   :members:

Miscellaneous functions
^^^^^^^^^^^^^^^^^^^^^^^

.. doxygengroup:: MiscRebFunctions
   :members:


