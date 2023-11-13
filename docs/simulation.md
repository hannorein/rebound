# Life cycle

These pages describe the C structure `reb_simulation` and the python class `rebound.Simulation`.
Because the python class is an abstraction of the C structure, we describe their common features in one place.
We will refer to both the C structure and the python object interchangeably as the *simulation structure* and *simulation object*.
The simulation structure contains all the configuration, status and particle data of one REBOUND simulation.
It's the one structure you will work with most when using REBOUND.



=== "C"
    The `reb_simulation_create()` function allocate memory for a `reb_simulation` structure and also initialize all variables to their default value.
    If you want to avoid a memory leak, you need to free the simulation when you no longer need it.
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... do work ... 
    reb_simulation_free(r);
    ```
    The call to `reb_simulation_free()` frees all data associated with the simulation such as particles.
    It will also free the memory for the simulation structure itself.

=== "Python"
    When you create a new object of the class `rebound.Simulation`, REBOUND will allocate memory for the object and also initialize all variables to their default value.

    Python automatically releases all the memory after the last reference to the object is gone:
    ```python
    sim = rebound.Simulation()
    # ... do work ...
    sim = None  # This will allow python to free the memory
    ```
    In general, you do not need to do this manually. 
    Python will loose the last reference to the simulation object at the end of the current variable scope (e.g. function). 

    !!! Danger
        The following code keeps a pointer to the `particles` array after the last reference to the simulation is gone.
        Because the memory associated with the `particles` array is freed when the simulation is freed, this will lead to a segmentation fault.

        ```python
        sim = rebound.Simulation()
        sim.add(m=1)
        particles = sim.particles
        sim = None           # free simulation
        print(sim.particles) # segmentation fault
        ```

There are several instances where you need to initialize a simulation's root boxes:

- If you use a tree code for collision detection or for calculating gravity.
- If you want to use open, periodic or shear-periodic boundary conditions.

Initializing root boxes is done after the simulation is created:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    double size = 100.;
    reb_simulation_configure_box(r, size, 1, 2, 3);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    size = 100.
    sim.configure_box(size, 1, 2, 3);
    ```
In the above example, there is one root box in the x direction, there are two in the y direction, and three in the z direction. 
In most cases you want exactly one root box in each direction.

