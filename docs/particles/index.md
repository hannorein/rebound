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

`#!c uint32_t hash`
:   integer or hash value used to identify the particle
    
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

Alternatively you can assign a hash value to particles and access them using the following syntax: 
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m", 1.); 
    r->particles[0].hash = reb_hash("star");
    reb_simulation_add_fmt(r, "a", 1.); 
    r->particles[1].hash = reb_hash("planet1");
    struct reb_particle* p = reb_simulation_particle_by_hash(r, reb_hash("planet1"));
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1., hash="star")
    sim.add(a=1., hash="planet1")
    p = sim.particles["planet1"]
    ```

