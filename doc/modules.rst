Available modules
=================

REBOUND is very modular. You have the choice between different gravity, collision, boundary and integrator modules. It is also possible to implement completely new modules with minimal effort. In the new version of REBOUND, modules are chosen at runtime by setting flags in the `reb_simulation` structure. 

The following sections list the available modules that come with REBOUND.

Gravity solvers
---------------
 
=======================  ============================================ 
Module name               Description
=======================  ============================================ 
REB_GRAVITY_COMPENSATED   Direct summation with compensated summation, O(N^2), default
REB_GRAVITY_NONE          No self-gravity
REB_GRAVITY_BASIC         Direct summation, O(N^2)
REB_GRAVITY_TREE          Oct tree, Barnes & Hut 1986, O(N log(N))
REB_GRAVITY_OPENCL        (upgrade to REBOUND 2.0 still in progress) Direct summation, O(N^2), but accelerated using the OpenCL framework.
REB_GRAVITY_FFT           (upgrade to REBOUND 2.0 still in progress) Two dimensional gravity solver using FFTW, works in a periodic box and the shearing sheet. 
=======================  ============================================ 


Collision detection algorithms
----------------------------

=======================  ============================================ 
Module name               Description
=======================  ============================================ 
REB_COLLISION_NONE        No collision detection, default
REB_COLLISION_DIRECT      Brute force collision search, O(N^2), checks for instantaneous overlaps only 
REB_COLLISION_LINE        Brute force collision search, O(N^2), checks for overlaps that occured during the last timestep assuming particles travelled along straight lines
REB_COLLISION_TREE        Oct tree, O(N log(N))
REB_COLLISION_SWEPP       (upgrade to REBOUND 2.0 still in progress) Plane sweep algorithm, ideal for low dimensional  problems, O(N) or O(N^1.5) depending on geometry 
=======================  ============================================ 


Boundary conditions
-------------------

=======================  ============================================ 
Module name               Description
=======================  ============================================ 
REB_BOUNDARY_NONE         Dummy. Particles are not affected by boundary conditions, default
REB_BOUNDARY_OPEN         Particles are removed from the simulation if they leaves the box.
REB_BOUNDARY_PERIODIC     Periodic boundary conditions. Particles are reinserted on the other side if they cross the box boundaries. You can use an arbitrary number of ghost-boxes with this module.
REB_BOUNDARY_SHEAR        Shear periodic boundary conditions. Similar to periodic boundary conditions, but ghost-boxes are moving with constant speed, set by the shear.
=======================  ============================================ 
 

Integrators
-----------

==========================  ============================================ 
Module name                 Description
==========================  ============================================ 
REB_INTEGRATOR_IAS15        IAS15 stands for Integrator with Adaptive Step-size control, 15th order. It is a vey high order, non-symplectic integrator which can handle arbitrary (velocity dependent) forces and is in most cases accurate down to machine precision. IAS15 can integrate variational equations. Rein & Spiegel 2015, Everhart 1985, default
REB_INTEGRATOR_WHFAST       WHFast is the integrator described in Rein & Tamayo 2015, it's a second order symplectic Wisdom Holman integrator with 11th order symplectic correctors. It is extremely fast and accurate, uses Gauss f and g functions to solve the Kepler motion and can integrate variational equations. The user can choose between Jacobi, Democratic Heliocentric, and WHDS (Hernandez and Dehnen, 2017) coordinates. 
REB_INTEGRATOR_JANUS        Janus is a bit-wise time-reversible high-order symplectic integrator using a mix of floating point and integer arithmetic. This integrator is still in an experimental stage and will be discussed in an upcoming paper. 
REB_INTEGRATOR_EULER        Euler scheme, first order
REB_INTEGRATOR_LEAPFROG     Leap frog, second order, symplectic
REB_INTEGRATOR_SEI          Symplectic Epicycle Integrator (SEI), mixed variable symplectic integrator for the shearing sheet, second order, Rein & Tremaine 2011
REB_INTEGRATOR_MERCURIUS    A hybrid integrator very similar to the one found in MERCURY. It uses WHFast for long term integrations but switches over smoothly to IAS15 for close encounters. This is a new integrator and might contain bugs. 
REB_INTEGRATOR_HERMES       A hybrid symplectic integrator that uses WHFast for long term integrations but switches over to IAS15 for close encounters. Note that this integrator is new and might contain bugs. In most cases the MERCURIUS integrator is perfomring better.
==========================  ============================================ 


