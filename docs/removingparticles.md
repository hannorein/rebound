# Removing particles

## Removing all particles

You can remove all particles with the following code:

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... add particles ...
    reb_simulation_remove_all_particles(r);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    del sim.particles
    ```
## Remove particle by index
Each particle in a REBOUND simulation can be uniquely identified with its position in the `particles` array, it's **index**.
You can remove a particle using this index as shown in the following code:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m", 1.); // star, index=0
    reb_simulation_add_fmt(r, "a", 1.); // planet 1, index=1
    reb_simulation_add_fmt(r, "a", 2.); // planet 2, index=2
    reb_simulation_remove_particle(r, 1, 1); // removes planet 1 (index 1)
    ```
    The first argument of `reb_simulation_remove_particle` is the simulation from which you want to remove the particle.
    The second argument is the index of the particle.
    The third argument determines if you want to keep the particle array sorted.
    In most cases you want to (set the argument to 1). 
    For simulation with many particles (millions), this might be slow. In that case set this argument to 0.

    The function returns 1 if the particle was successfully removed, and 0 if the index was out of range.

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1.) // star, index=0
    sim.add(a=1.) // planet 1, index=1
    sim.add(a=2.) // planet 2, index=2
    sim.remove(1)
    ```
    The `remove` function accepts an optional argument `keep_sorted`. 
    It determines if you want to keep the particle array sorted.
    In most cases you want to (set the argument to `True`, the default). 
    For simulation with many particles (millions), this might be slow. In that case set this argument to `False`.

## Remove particle by hash
In addition to the index, you can also identify particles by hash.
This is useful when you want to make sure you can uniquely identify particles in simulations where you constantly add or remove particles. 
If a particle has a hash, you can remove it as shown here:

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m", 1.); 
    r->particles[0].hash = reb_hash("star");
    reb_simulation_add_fmt(r, "a", 1.); 
    r->particles[1].hash = reb_hash("planet1");
    reb_simulation_add_fmt(r, "a", 2.); 
    r->particles[2].hash = reb_hash("planet2");
    reb_simulation_remove_particle_by_hash(r, reb_hash("planet1"), 1);
    ```
    The syntax of the function is the same as for `reb_simulation_remove_particle` except you pass the hash instead of the index.

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1., hash="star")
    sim.add(a=1., hash="planet1")
    sim.add(a=2., hash="planet2")
    sim.remove(hash="planet1")
    ```

