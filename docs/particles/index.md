# Particle structure
A particle is represented by the `reb_particle` structure in C.
The python class `Particle` is an abstraction of the `reb_particle` structure in C.
We will refer to both the C structure and the python object interchangeably as the *particle structure* and *particle object*.

The particle object contains the following variables which can be directly manipulated:

`#!c double m`
:   mass

`#!c double r`
:   physical radius of the particle

`#!c double x, y, z, vx, vy, vz`
:   position and velocity coordinates

`#!c const char* name`
:   The name of the particle as a NULL terminated string. If the particle is part of a REBOUND simulation, then the memory of the name is managed by REBOUND. Do not set the name of a particle directly if it is in a REBOUND simulation. Use `reb_particle_set_name()` instead.
    
You can create a particle object which is not part of a REBOUND simulation.
=== "C"
    ```c
    struct reb_particle p = {.m=1., x=0., vy=0.};
    ```

=== "Python"
    ```python
    p = rebound.Particle(m=1., x=0., vy=0.)
    ```

However, in most cases you will work with particles which have been added to a REBOUND simulation.
You then access the particle using the simulation's `particles` array:

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation, add particles ...
    r->particles[0].x = 1;
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... setup simulation, add particles ... 
    sim.particles[0].x = 1
    ```

Alternatively you can assign a name to particles and access them using the following syntax: 
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m name", 1., "star"); 
    reb_simulation_add_fmt(r, "a name", 1., "planet1"); 
    struct reb_particle* p = reb_simulation_get_particle_by_name(r, "planet1");
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1., name="star")
    sim.add(a=1., name="planet1")
    p = sim.particles["planet1"]
    ```

