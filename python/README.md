IAS15 - Python bindings
===========================================================================

Contributors
------------
* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>

Workflow
------------
This folder contains python bindings for the high order, high accuracte IAS15 integrator. 
It does not support all features of REBOUND. These bindings are for those users who would
like to use python instead of C. Since the backend of the integrator is running in C,
one can expect C-like performance.

The work flow is very simple, everything more complicted is left up to you an external package.

1) Setup the gravitational constant G and the initial timestep dt.
2) Setup particles in Cartesian coordinates.
3) Pass the particles to the REBOUND module
4) Integrate for a fixed number of steps or for a amount fixed time.
5) Read back the particle positions in Cartesian coordinates.

Installation
------------
The python module makes use of a shared dynamic library and ctypes. You only need to 
compile the shared library (this will generate libias15.so):

    make

Then, look at the examples in the differen folders. All examples use the rebound module 
(rebound.py) which encapsulates the interface to the code.
