# Visualization

Starting with version 4, REBOUND includes new real-time interactive 3D visualizations. 
The code for these visualization is built-upon the previous OpenGL visualizations that came with the C version of REBOUND. 
However, the new visualization feature has several advantages

- There are zero dependecies. No need to install GLFW or other libraries.
- The visualizations work on Linux, MacOS, Windows, including mobile devices.
- The visualizations work both for the C and python version of REBOUND.
- You can use visualizations for simulation running on remote servers.

This page describes how to use these visualizations and the technology that makes this possible behind the scenes.

## Using G=1

As an example, suppose we use $G=1$, have a central object of mass $M$, and a test particle orbiting on a circular orbit at a distance $a=1$. 
This scenario can be setup with the following code:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt("m", 1.);
    reb_simulation_add_fmt("a", 1.);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(a=1.)
    ```
