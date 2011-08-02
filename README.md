REBOUND - A collisional N-body code
========================================
2011 ISIMA project

Contributors
------------
* Hanno Rein, Institute for Advanced Study (IAS), Princeton, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>
  
Available modules
-----------------
* Gravity
   - No self-gravity
   - Direct summation, O(N^2)
   - Oct tree, Barnes & Hut 1986, O(N log(N))
* Integrators
   - Euler, first order
   - Leap frog, second order, symplectic
   - Mixed variable symplectic integrator for 1/r potential, second order, Wisdom & Holman 1991, Kinoshita et al 1991
   - Symplectic Eplicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011
* Collision detection
   - No collision detection
   - Direct nearest neighbour search, O(N^2)
   - Oct tree, O(N log(N))
* Output/Visualization
   - Standard ASCII or binary output 
   - Real-time, 3D OpenGL visualization

Other features
--------------
* The code is written entirely in C. It is standard compliant to C99.
* Parallilized with OpenMP (for shared memory systems)
* Parallilized with MPI using an essentrial tree for gravity and collisions (for distributed memory systems)
* No libraries are needed. The use of OpenGL/GLUT/libpng for visualization is optional.
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/rebound.
* No configuration is needed to run any of the example problems. Just type 'make && ./nbody' to run them.
* Different modules are easily interchangable by one line in the Makefile.
  

How to compile and run
----------------------

**For the impatient**

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./nbody

**Long description**

rebound is very easy to use. To get started, download the latest version of the code from github. If you are familiar with `git`, you can clone the project and keep up-to-date with the latest developments. Otherwise, you can also simply download a snapshot of the repository. 

In the main directory, you find a sub-directory called `src` which contains the source code and a directory called `examples` with various example problems. To compile one of the example, go have to go to that direcotry, for example:

    cd examples/shearing_sheet/

Then, type

    make

This will do the following things    

* It sets various environment variables. These determine the compiler optimization flags and which libraries should be included. In the example `shearing_sheet` the variables are:
   - `OPT=-O3`. This sets the additional compiler flag `-O3` which optimizes the code for speed.
   - `QUADRUPOLE=0`. This disables the calculation of quadrupole moments in the tree code. The simulation is faster, but less accurate.
   - `OPENGL=1`. This enables real-time OpenGL visualizations and requires both OpenGL and GLUT libraries to be installed. This should work without any further adjustments on any Mac that has Xcode installed. On Linux both libraries must be installed in `/usr/local/`. You can change the default search paths for libraries in the file `src/Makefile`. 
   - `MPI=0`. This disables parallelization with MPI.
   - `OPENMP=1`. This enables parallelization with OpenMP. The number of threads can be set with an environment variable at runtime, e.g.: `export OMP_NUM_THREADS=8`.
* It creates symbolic links to the active modules. This allows you to choose from different gravity solvers, boundary conditions, integrators and collision solvers. For example, to change the gravity solver from using a tree to direct summation you could change `gravity_tree.c` to `gravity_direct.c`. 
* It creates a symbolic link to the current problem file. A problem files contains the initial conditions and the output routines for the current problem. You do not need to change any file in `src/` to create a new problem unless you want to do something very special. This keeps the inital conditions and the code itself cleanly separated.
* It compiles the code and copies the binary into the current directory.

You can also create a documentation with `doxygen` based on the current choice of modules by typing `make doc`. This requires `doxygen` to be installed. The documentation will be generated in the directory `doc/html/`.

To run the code, simply type

    ./nbody

If you want to create your own problem, just copy one of the example directories and modify `problem.c` and `Makefile` accordingly.     


License
-------
rebound is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

rebound is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rebound.  If not, see <http://www.gnu.org/licenses/>.

Acknowledgements
----------------
When you use this code or parts of this code for results presented in a scientific publication, we would greatly appriciate a citation to Rein and Liu (in preparation) and an acknoledgement of the form: 

_Calculations in this paper made use of the collisional N-body code rebound which can be downloaded freely at http://github.com/hannorein/rebound._

Also, please send us a copy of your paper so that we can keep track of all publications that made use of the code.
