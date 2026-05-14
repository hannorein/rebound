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
