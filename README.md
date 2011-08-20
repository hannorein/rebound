REBOUND - An open-source multi-purpose N-body code for collisional dynamics
===========================================================================

Contributors
------------
* Hanno Rein, Institute for Advanced Study (IAS), Princeton, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>


Screenshot
---------- 
 
![Tree structure in REBOUND](/hannorein/rebound/raw/master/doc/images/screenshot_shearingsheet.png) 

Available modules
-----------------

REBOUND is extremely modular. You have the choice between different gravity, collision, boundary and integration modules. It is also possible to implement completely new modules with minimal effort. Modules are chosen by setting symbolic links. Thus, there is no need to run a configure script. For example, there is one link `gravity.c` that points to one of the gravity modules `gravity_*.c`. The symbolic links are set in the problem makefile (see below).

This setup allows you to work on multiple projects at the same time using different modules. When switching to another problem, nothing has to be set-up and the problem can by compiled by simply typing `make` in the corresponding directory (see below).

### Gravity ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>gravity_none.c       </pre></td>
     <td>No self-gravity</td></tr>
  <tr><td><pre>gravity_direct.c     </pre></td>
     <td>Direct summation, O(N^2)</td></tr>
  <tr><td><pre>gravity_tree.c       </pre></td>
     <td>Oct tree, Barnes & Hut 1986, O(N log(N))</td></tr>
</table>

### Collision detection ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>collisions_none.c    </pre></td>
     <td>No collision detection</td></tr>
  <tr><td><pre>collisions_direct.c  </pre></td>
     <td>Direct nearest neighbor search, O(N^2)</td></tr>
  <tr><td><pre>collisions_tree.c    </pre></td>
     <td>Oct tree, O(N log(N))</td></tr>
  <tr><td><pre>collisions_sweep.c   </pre></td>
     <td>Plane sweep algorithm, ideal for low dimensional problems, O(N) or O(N^1.5) depending on geometry</td></tr>
</table>

### Integrators ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>integrator_euler.c   </pre></td>
     <td>Euler scheme, first order</td></tr>
  <tr><td><pre>integrator_leapfrog.c</pre></td>
     <td>Leap frog, second order, symplectic</td></tr>
  <tr><td><pre>integrator_wh.c      </pre></td>
     <td>Wisdom-Holman Mapping, mixed variable symplectic integrator for the Kepler potential, second order, Wisdom & Holman 1991, Kinoshita et al 1991</td></tr>
  <tr><td><pre>integrator_sei.c     </pre></td>
     <td>Symplectic Epicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011</td></tr>
</table>

### Boundaries ###
<table>
  <tr><th>Module name</th>
     <th>Description</th></tr>
  <tr><td><pre>boundaries_open.c    </pre></td>
     <td>Particles are removed from the simulation if they leaves the box.</td></tr>
  <tr><td><pre>boundaries_periodic.c</pre></td>
     <td>Periodic boundary conditions. Particles are reinserted on the other side if they cross the box boundaries. You can use an arbitrary number of ghostboxes with this module.</td></tr>
  <tr><td><pre>boundaries_shear.c   </pre></td>
     <td>Shear periodic bounary conditions. Similar to periodic boundary conditions, but ghostboxes are moving with constant speed, set by the shear.</td></tr>
</table>


Other features worth mentioning
-------------------------------
* Real-time, 3D OpenGL visualization.
* The code is written entirely in C. It conforms to the ISO standard C99.
* Parallelized with OpenMP (for shared memory systems).
* Parallelized with MPI using an essential tree for gravity and collisions (for distributed memory systems).
* No libraries are needed. The use of OpenGL/GLUT/libpng for visualization is optional.
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/rebound.
* No configuration is needed to run any of the example problems. Just type `make && ./nbody` in the problem directory to run them.
* Standard ASCII or binary output routines. 
* Different modules are easily interchangeable by one line in the Makefile.
  

How to download, compile and run REBOUND
----------------------------------------

### For the impatient ###
Simply copy and paste this line to your terminal and press enter

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./nbody

or if you do not have git installed

    wget https://github.com/hannorein/rebound/tarball/master -O- | tar xvz && cd hannorein-rebound-*/examples/shearing_sheet/ && make && ./nbody

### For the patient ###
REBOUND is very easy to install and use. To get started, download the latest version of the code from github. If you are familiar with `git`, you can clone the project and keep up-to-date with the latest developments. Otherwise, you can also simply download a snapshot of the repository as a tar or zip file at http://github.com/hannorein/rebound. There is a download bottom at the top right. 

In the main directory, you find a sub-directory called `src` which contains the bulk parts of the  source code and a directory called `examples` with various example problems. To compile one of the example, you have to go to that directory, for example:

    cd examples/shearing_sheet/

Then, type

    make

This will do the following things    

* It sets various environment variables. These determine settings like the compiler optimization flags and which libraries are included (see below). 
* It creates symbolic links to the active modules. This allows you to choose from different gravity solvers, boundary conditions, integrators and collision solvers. For example, to change the gravity solver from using a tree to direct summation you could change `gravity_tree.c` to `gravity_direct.c`. 
* It creates a symbolic link to the current problem file. Each problem file contains the initial conditions and the output routines for the current problem. You do not need to change any file in `src/` to create a new problem unless you want to do something very special. This keeps the initial conditions and the code itself cleanly separated.
* It compiles the code and copies the binary into the current directory.

If something goes wrong, it is most likely the visualization module. You can turn it off by deleting the line which contains `OPENGL` in the makefile. Of course, you will not see much unless you put in some extra work to visualize the results.

You can also create a documentation based on the current choice of modules by typing `make doc`. However, this requires the documentation generator `doxygen` to be installed. The documentation will be generated in the directory `doc/html/`.

To finally run the code, simply type

    ./nbody

A window should open and you will see a simulation running in real time. The setup simulates the rings of Saturn and uses a local shearing sheet approximation. Have a look at the other examples too and you will quickly get an impression of what REBOUND can do. 

If you want to create your own problem, just copy one of the example directories or the template in the `problems` directory. Then simply modify `problem.c` and `Makefile` accordingly.  

### Environment variables ###
The makefile in each problem directory sets various environment variables. These determine the compiler optimization flags, the libraries included and basic code settings. Let us look at one of the examples `shearing_sheet` in more detail. 

- `export QUADRUPOLE=0`. This disables the calculation of quadrupole moments for each cell in the tree. The simulation is faster, but less accurate.
- `export OPENGL=1`. This enables real-time OpenGL visualizations and requires both OpenGL and GLUT libraries to be installed. This should work without any further adjustments on any Mac which has Xcode installed. On Linux both libraries must be installed in `/usr/local/`. You can change the default search paths for libraries in the file `src/Makefile`. 
- `export MPI=0`. This disables parallelization with MPI.
- `export OPENMP=1`. This enables parallelization with OpenMP. The number of threads can be set with an environment variable at runtime, e.g.: `export OMP_NUM_THREADS=8`.
- `export CC=icc`. This flag can be used to override the default compiler. The default compilers are `gcc` for the sequential and `mpicc` for the parallel version. 
- `export LIB=`. Additional search paths for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 
- `export OPT=-O3`. This sets the additional compiler flag `-O3` and optimizes the code for speed. Additional search paths to header files for external libraries (such as OpenGL, GLUT and LIBPNG) can be set up using this variable. 

All of these variables are read by the main makefile in the `src/` directory. The `OPENGL` variable, for example, is used to determine if the OpenGL and GLUT libraries should be included. If the variable is `1` the makefile also sets a pre-compiler macro with `-DOPENGL`. Note that because OPENGL is incompatible with MPI, when MPI is turned on (set to 1), OPENGL is automatically turned off (set to 0) in the main makefile.


### User-defined functions in the problem.c file ###
The problem.c file contains at least four functions. You do not need to implement all of them, a dummy might be enough. 

#### void problem_init(int argc, char* argv[]) ####
This routine is where you read command line arguments and set up your initial conditions. REBOUND does not come with a built-in functionality to read configuration files at run-time. You should see this as a feature. In REBOUND, you have one `problem.c` file for each problem. Thus, everything can be set within this file. There are, of course, situation in which you want to do something like a parameter space survey. In almost all cases, you vary a few parameters but rarely more, say 5. You can easily read these parameters from the command line.
 
Here is one example that reads in the first argument given to rebound as the box-size and sets a default value when no value is given:

```c
if (argc>1){
	boxsize = atof(argv[1]);
}else{
	boxsize = 100;
}
```

If you are still convinced that you need a configuration file, you are welcome to implement it yourself. This function is where you want to do that.    

#### void problem_inloop() ####
This function is called once per time-step. It is called at the end of the K part of the DKD time-stepping scheme. This is where you can implement all kind of things such as additional forces onto particles. 

The following lines of code, for example, implement the Poynting Robertson drag force on each particle except the first one (which is the star in this example):

```c
double alpha = 1e-4;
for (int i=1;i<N;i++){
	double x = particles[i].x;
	double y = particles[i].y;
	double z = particles[i].z;
	double r2 = (x*x + y*y + z*z); 
	particles[i].vx -= dt * particles[i].vx*alpha/r2;
	particles[i].vy -= dt * particles[i].vy*alpha/r2;
	particles[i].vz -= dt * particles[i].vz*alpha/r2;
}
```

#### void problem_output() ####
This function is called at the beginning of the simulation and at the end of each time-step. You can implement your output routines here. Many basic output functions are already implemented in REBOUND. See `output.h` for more details. The function `output_check(odt)` can be used to easily check if an output is needed after a regular interval. For example, the following code snippet outputs some timing statistics to the console every 10 time-steps:

```c
if (output_check(10.*dt)){
	output_timing();
}
```    
 
#### void problem_finish() ####
This function is called at the end of the simulation, when t >= tmax. This is the last chance to output any quantities before the program ends.


Support and contributions
-------------------------
We offer limited support for REBOUND. If you encounter any problems, just send us an e-mail with as much details as possible and include your problem.c and makefile. Please make sure you are using the latest version of REBOUND that is available on github. 

REBOUND is open source and you are strongly encouraged to contribute to this project if you are using it. Please contact us and we will give you permission to push directly to the public repository. 


License
-------
REBOUND is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

REBOUND is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with REBOUND.  If not, see <http://www.gnu.org/licenses/>.


Acknowledgments
---------------
When you use this code or parts of this code for results presented in a scientific publication, we would greatly appreciate a citation to Rein and Liu (in preparation) and an acknowledgment of the form: 

_Simulations in this paper made use of the collisional N-body code REBOUND which can be downloaded freely at http://github.com/hannorein/rebound._

Also, please send us a copy of your paper so that we can keep track of all publications that made use of the code.
