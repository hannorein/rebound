2011 ISIMA project on N-body simulations
========================================

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
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/nbody.
* No configuration is needed to run any of the example problems. Just type 'make && ./nbody' to run them.
* Different modules are easily interchangable by one line in the Makefile.
  
License
-------
nbody is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

nbody is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with nbody.  If not, see <http://www.gnu.org/licenses/>.

Acknowledgements
----------------
When you use this code or parts of this code for results presented in a scientific publication, we would greatly appriciate a citation to #TBD# and an acknoledgement of the form: 

_Calculations in this paper made use of the collisional N-body code #TBD# which can be downloaded freely at http://github.com/hannorein/nbody._

Also, please send us a copy of your paper so that we can keep track of all publications that made use of the code.
