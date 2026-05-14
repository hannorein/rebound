# Saving simulations to disk
You can use binary files to save simulations to a file and then later restore them from this file.
All the particle data and the current simulation states are saved. 
Below is an example on how to work with binary files.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation ...
    reb_simulation_integrate(r, 10); // integrate 
    reb_simulation_save_to_file(r, "snapshot.bin");
    reb_simulation_free(r); 

    struct reb_simulation* r2 = reb_simulation_create_from_file("snapshot.bin", 0);
    reb_simulation_integrate(r2, 20); // continue integration
    reb_simulation_free(r2); 
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    // ... setup simulation ...
    sim.integrate(10)
    sim.save_to_file("snapshot.bin")
    sim = None # Remove reference, allow python to release memory

    sim2 = rebound.Simulation("snapshot.bin")
    sim2.integrate(2) # continue integration
    sim2 = None 
    ```

Rather than using one file for one snapshot of a simulation, you can also use a [Simulationarchive](simulationarchive.md).
A Simulationarchive is a collection of simulation snapshots stored in one binary file. 



