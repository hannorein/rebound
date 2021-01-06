Welcome to REBOUND!
===================
.. image:: images/dense.png
   :width: 362px
   :height: 350px
   :align: right
.. image:: images/disc.png
   :width: 362px
   :height: 350px
   :align: right
.. image:: images/outersolarsystem.png
   :width: 362px
   :height: 350px
   :align: right
.. image:: images/shearingsheet.png
   :width: 362px
   :height: 350px
   :align: right

REBOUND is an N-body integrator, i.e. a software package that can integrate the motion of particles under the influence of gravity. The particles can represent stars, planets, moons, ring or dust particles. REBOUND is very flexible and can be customized to accurately and efficiently solve many problems in astrophysics.  An incomplete feature list of REBOUND:

* Symplectic integrators (WHFast, WHFastHelio, SEI, LEAPFROG)
* High order symplectic integrators (SABA, WH Kernel methods)
* High accuracy non-symplectic integrator with adaptive timestepping (IAS15)
* Support for collisional/granular dynamics, various collision detection routines
* The code is written entirely in C, conforms to the ISO standard C99 and can be used as a thread-safe shared library
* Easy-to-use Python module, installation in 3 words: `pip install rebound`
* Extensive set of example problems in both C and Python
* Real-time, 3D OpenGL visualization (C version)
* Parallelized with OpenMP (for shared memory systems)
* Parallelized with MPI using an essential tree for gravity and collisions (for distributed memory systems)
* No libraries are needed, use of OpenGL/glfw3 for visualization is optional
* The code is fully open-source and can be downloaded freely from http://github.com/hannorein/rebound
* No configuration is needed to run any of the example problems. Just type `make && ./rebound` in the problem directory to run them
* Comes with standard ASCII or binary output routines 
* Different modules are easily interchangeable at runtime



How to use REBOUND - a quick introduction
-----------------------------------------
You can call REBOUND from C or Python. Which programming language you want to use depends on your taste and your specific application. In short: If you simply want to setup a few particles such as a planetary system, visualize it with a WebGL widget, and integrate it with the high order integrator IAS15 or the new symplectic integrator WHFast then use the Python version. If you want to run large simulations with millions of particles, use an exotic integrator, use fast OpenGL visualizations, or make use of the distributed tree code then use the C version. 

All the computationally expensive parts of REBOUND are written in C. So even if you use the Python version, you'll end up with a very fast code.

To install the *Python version*, simply type the following command into a terminal::

    pip install rebound

To learn more about how to use REBOUND with Python have a look at the iPython/Jupyter tutorials at https://github.com/hannorein/rebound/blob/master/ipython_examples/

To install the *C version*, clone this repository, e.g. by simply copy-and-pasting the following command into your terminal::
    
    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

To learn more about how to use REBOUND with C, study the examples in the `examples/` directory and continue reading this file. You might also want to have a look at the `rebound.h` file in the `src/` directory which contains the API specifications. Last but not least, REBOUND is open source. If you want to know how something works, you can just look at the source code. And of course, you are welcome to e-mail any of the contributors with questions. We'll do our best to answer them quickly.

Note:  If you want to run REBOUND on Windows, the best way is likely to install the Windows Subsystem for Linux. After installing the gcc compiler, e.g., sudo apt-get install gcc, you should be able to install REBOUND and any python libraries by following the Linux/Mac installation instructions in this documentation. Unfortunately we do not have Windows installations ourselves, so we cannot actively support installation problems. Thanks to Keto /Zhang for finding this workaround.

YouTube Videos
--------------
Short YouTube videos describing various aspects of REBOUND are available at https://www.youtube.com/channel/UC2wonKI0wWwGi5-JqJtMsYQ/videos.

Contributors
------------
* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>
* Akihiko Fujii, National Astronomical Observatory of Japan/University of Tokyo, Tokyo, <akihiko.fujii@nao.ac.jp>
* Dan Tamayo, University of Toronto, <dtamayo@cita.utoronto.ca>
* Ari Silburt, Penn State University <ajs725@psu.edu>
* and many others!

REBOUND is open source. You are invited to contribute to this project if you are using it. Please contact any of the authors above if you have any questions.


Papers
------
There are several papers describing the functionality of REBOUND. 

1. Rein & Liu 2012 (Astronomy and Astrophysics, Volume 537, A128) describes the code structure and the main feature including the gravity and collision routines for many particle systems. http://adsabs.harvard.edu/abs/2012A%26A...537A.128R 

2. Rein & Tremaine 2011 (Monthly Notices of the Royal Astronomical Society, Volume 415, Issue 4, pp. 3168-3176) describes the Symplectic Epicycle integrator for shearing sheet simulations. https://ui.adsabs.harvard.edu/abs/2011MNRAS.415.3168R 

3. Rein & Spiegel 2015 (Monthly Notices of the Royal Astronomical Society, Volume 446, Issue 2, p.1424-1437) describes the versatile high order integrator IAS15 which is now part of REBOUND. http://adsabs.harvard.edu/abs/2015MNRAS.446.1424R

4. Rein & Tamayo 2015 (Monthly Notices of the Royal Astronomical Society, Volume 452, Issue 1, p.376-388) describes WHFast, the fast and unbiased implementation of a symplectic Wisdom-Holman integrator for long term gravitational simulations. http://adsabs.harvard.edu/abs/2015MNRAS.452..376R

5. Rein & Tamayo 2016 (Monthly Notices of the Royal Astronomical Society, Volume 459, Issue 3, p.2275-2285) develop the framework for second order variational equations. https://ui.adsabs.harvard.edu/abs/2016MNRAS.459.2275R

6. Rein & Tamayo 2017 (Monthly Notices of the Royal Astronomical Society, Volume 467, Issue 2, p.2377-2383) describes the Simulation Archive for exact reproducibility of N-body simulations. https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.2377R

7. Rein & Tamayo 2018 (Monthly Notices of the Royal Astronomical Society, Volume 473, Issue 3, p.3351–3357) describes the integer based JANUS integrator. https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.3351R

8. Rein, Hernandez, Tamayo, Brown, Eckels, Holmes, Lau, Leblanc & Silburt 2019 (Monthly Notices of the Royal Astronomical Society, Volume 485, Issue 4, p.5490-5497) describes the hyrbid symplectic integrator MERCURIUS. https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.5490R

9. Rein, Tamayo & Brown 2019 (Monthly Notices of the Royal Astronomical Society, Volume 489, Issue 4, November 2019, Pages 4632-4640) describes the implementation of the high order symplectic intergators SABA, SABAC, SABACL, WHCKL, WHCKM, and WHCKC. https://ui.adsabs.harvard.edu/abs/2019MNRAS.489.4632R


Acknowledgments
---------------
If you use this code or parts of this code for results presented in a scientific publication, please send us a copy of your paper so that we can keep track of all publications that made use of the code. We would greatly appreciate a citation as well. The simplest way to find the citations relevant to the specific setup of your REBOUND simulation is: 


    sim = rebound.Simulation()
    -your setup-
    sim.cite()


License
-------
REBOUND is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

REBOUND is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with REBOUND.  If not, see <http://www.gnu.org/licenses/>.


Table of Contents
-----------------
.. toctree::
   :maxdepth: 2

   
   self
   changelog
   modules
   quickstart
   examples
   c_api
   python_api

