# Integrators

![type:video](https://www.youtube.com/embed/QW5a-iH62dQ)

Numerical integrators are the backbone of any N-body package. 
A numerical integrator evolves particles forward in time, one timestep at a time.
To do that, the integrator needs to know the current position and velocity coordinates of the particles, and the equations of motion which come in the form of a set of ordinary differential equations.

Because an exact solution to these differential equations is in general unknown, each integrator attempts to approximate the true solution numerically. 
Different integrators do this differently and each of them has some advantages and some disadvantages. 
Each of the built-in integrators of REBOUND is described in this section.

## Custom Integrator
Starting with version 4.6.1, REBOUND has an API which allows users to easily add their own integrator without changing REBOUND itself. 
This is particularly helpful for developing new numerical methods: you have complete control over the integration step, but can fall back to many features of REBOUND such as setting up orbits, calculating orbital parameters, calculating energy errors, etc. 
Furthermore, if your custom integrator needs to store any data, this can be done automatically with Simulationarchive snapshpts, allowing you to easily restart simulations. 

To enable a custom integrator, you must set `integrator` and setup various function pointers in the `ri_custom` struct. 
This feature is intended to be used from C and we do currently not provide a python interface.
There is one [basic C example](../c_examples/custom_integrator) and one more [advanced C example](../c_examples/custom_integrator_with_data) which show the usage.


