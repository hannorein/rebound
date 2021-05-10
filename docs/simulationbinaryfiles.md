# Saving simulations to disk
You can use binary files to save simulations to a file and then later restore them from this file.
All the particle data and the current simulation states are saved. 
Below is an example on how to work with binary files.

=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... setup simulation ...
    reb_integate(r, 10); // integrate 
    reb_output_binary(r, "snapshot.bin");
    reb_free_simulation(r); 

    struct reb_simulation* r2 = reb_create_simulation_from_binary("snapshot.bin);
    reb_integate(r2, 20); // continue integration
    reb_free_simulation(r2); 
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    // ... setup simulation ...
    sim.integate(10)
    sim.save("snapshot.bin")
    sim = None # Remove reference, allow python to release memory

    sim2 = rebound.Simulation("snapshot.bin")
    sim2.integrate(2) # continue integration
    sim2 = None 
    ```

Rather than using one file for one snapshot of a simulation, you can also use a [Simulation Archive](simulationarchive.md).
A Simulation Archive is a collection of simulation snapshots stored in one binary file. 



