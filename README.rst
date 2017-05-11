REBOUND - An open-source multi-purpose N-body code
==================================================

.. image:: http://img.shields.io/badge/rebound-v3.4.0-green.svg?style=flat
    :target: http://rebound.readthedocs.org
.. image:: https://badge.fury.io/py/rebound.svg
    :target: https://badge.fury.io/py/rebound
.. image:: http://img.shields.io/badge/license-GPL-green.svg?style=flat 
    :target: https://github.com/hannorein/rebound/blob/master/LICENSE
.. image:: http://img.shields.io/travis/hannorein/rebound/master.svg?style=flat 
    :target: https://travis-ci.org/hannorein/rebound/
.. image:: https://coveralls.io/repos/hannorein/rebound/badge.svg?branch=master&service=github 
    :target: https://coveralls.io/github/hannorein/rebound?branch=master
.. image:: http://img.shields.io/badge/arXiv-1110.4876-green.svg?style=flat 
    :target: http://arxiv.org/abs/1110.4876
.. image:: http://img.shields.io/badge/arXiv-1409.4779-green.svg?style=flat 
    :target: http://arxiv.org/abs/1409.4779
.. image:: http://img.shields.io/badge/arXiv-1506.01084-green.svg?style=flat 
    :target: http://arxiv.org/abs/1506.01084
.. image:: http://img.shields.io/badge/arXiv-1603.03424-green.svg?style=flat 
    :target: http://arxiv.org/abs/1603.03424 
.. image:: https://readthedocs.org/projects/pip/badge/?version=latest
    :target: http://rebound.readthedocs.org/
.. image:: https://img.shields.io/badge/launch-binder-ff69b4.svg?style=flat
    :target: http://mybinder.org/repo/hannorein/rebound



FEATURES
--------

REBOUND is an N-body integrator, i.e. a software package that can integrate the motion of particles under the influence of gravity. The particles can represent stars, planets, moons, ring or dust particles. REBOUND is very flexible and can be customized to accurately and efficiently solve many problems in astrophysics.  An incomplete feature list of REBOUND:

* Symplectic integrators (WHFast, WHFastHelio, SEI, LEAPFROG)
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

One minute installation
-----------------------

You can install REBOUND with pip if you want to only use the python version of REBOUND::

    pip install rebound

Then, you can run a simple REBOUND simulation such as

.. code:: python

   import rebound
   sim = rebound.Simulation()
   sim.add(m=1.0)
   sim.add(m=1.0e-3, a=1.0)
   sim.integrate(1000.)
   sim.status()

If you want to use the C version of REBOUND simply copy and paste this line into your terminal (it won't do anything bad, we promise)::

    git clone http://github.com/hannorein/rebound && cd rebound/examples/shearing_sheet && make && ./rebound

 
Documentation
-------------
The full documentation with many examples, changelogs and tutorials can be found at

http://rebound.readthedocs.org

We're alway trying to improve REBOUND and extending the documention is high on our to-do list.
If you have trouble installing or using REBOUND, please open an issue on github and we'll try to help as much as we can.


Changelog
---------
For a changelog of the most important changes in recent updates, see https://github.com/hannorein/rebound/blob/master/changelog.rst 
