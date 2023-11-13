# Simulationarchive

The concepts behind the Simulationarchive are described in detail in [Rein & Tamayo 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.2377R/abstract).
Further examples of how to work with the Simulationarchive are provided in an [iPython](ipython_examples/Simulationarchive.ipynb) and [C example](c_examples/simulationarchive.md).

## Creating Simulationarchive snapshots

The following code shows how to manually append a Simulationarchive snapshot to a file.
If the file does not exist yet, the function outputs a new binary file. 
If the file already exists, the function will append a Simulationarchive snapshot to the existing file. 

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... work on simulation ...
    reb_simulation_save_to_file(r, "archive.bin");
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.save_to_file("archive.bin")
    ```
    You can pass the optional argument `delete_file=True` to delete the file if it already exists.
    By default, the function appends a snapshot to existing files.

Instead of manually outputting each snapshot, you can also automate this process as shown below.

### Regular time intervals
The following code automatically creates a Simulationarchive snapshot at regular intervals.
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... work on simulation ...
    reb_simulation_save_to_file_interval(r, "archive.bin", 10.);
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.save_to_file("archive.bin", interval=10.)
    ```

### Regular number of timesteps
The following code automatically creates a Simulationarchive snapshot after a fixed number of timesteps.
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... work on simulation ...
    reb_simulation_save_to_file_step(r, "archive.bin", 100);
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.save_to_file("archive.bin", step=100)
    ```
!!! Info
    This method is in general more reliable than the interval method.
    The reason is that the number of timesteps is an integer value whereas the time is a floating point number.
    If you run long simulations, you might encounter issues with finite floating point precision.
    This only affects very high accuracy simulation where you want to make sure the outputs occur exactly at the right timestep. 


### Regular wall-time intervals
The following code automatically creates a Simulationarchive snapshot after a fixed wall-time.
This is particularly useful for creating restart files when running long simulations.
The wall-time is given in seconds.
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... work on simulation ...
    reb_simulation_save_to_file_walltime(r, "archive.bin", 120.); // 2 minutes between snapshots
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.save_to_file("archive.bin", walltime=120) # 2 minutes
    ```

## Reading Simulationarchives

### Reading one snapshot
The following example shows how to read in a specific snapshot of a Simulationarchive.
If you pass a negative number for the snapshot, it will wrap around to the end of the Simulationarchive.
For example, the last snapshot in the file would have the index `-1`, the second to last `-2`, and so on.
=== "C"
    ```c
    struct reb_simulationarchive* archive = reb_simulationarchive_create_from_file("archive.bin");
    struct reb_simulation* r = reb_simulation_create_from_simulationarchive(archive, 12); // snapshot with index 12
    reb_simulationarchive_free(archive);
    // ... work on simulation ...
    ```

=== "Python"
    ```python
    sim = rebound.Simulation("archive.bin", snapshot=12)
    # ... work on simulation ...
    ```
    In python, `snapshot=-1` is the default. 
    Thus, `#!python sim = rebound.Simulation("archive.bin")` will create a new simulation from the last snapshot in the archive. 

