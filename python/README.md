IAS15 - Python bindings
===========================================================================

Contributors
------------
* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>

Workflow
------------
This folder contains python bindings for the high order, high accurate IAS15 integrator. 
It does not support all features of REBOUND. These bindings are for those users who would
like to use python instead of C. Since the backend of the integrator is running in C (the 
exact same code that REBOUND uses, one can expect C-like performance.

The work flow is very simple. 

1. Setup the gravitational constant G and the initial timestep dt.
2. Setup particles in Cartesian coordinates.
3. Pass the particles to the REBOUND module
4. Integrate for a fixed number of steps or for a amount fixed time.
5. Read back the particle positions in Cartesian coordinates.

Everything more complicated is left up to the user (e.g. input/output, calculating orbital elements, visualization). 

Installation
------------
The rebound python module makes use of a shared dynamic library and ctypes. You only need to 
compile the shared library once (this will generate libias15.so):

    make

Then, look at the examples in the subfolders in this directory. All examples should be 
self-explanatory. They all use the rebound module (rebound.py) which encapsulates the 
interface to the shared library. 

To run the examples, simply go into a directory and type:

    python problem.py
