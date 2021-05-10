# Simulation Archive

The concepts behind the SimulationArchive are described in detail in [Rein & Tamayo 2017](https://ui.adsabs.harvard.edu/abs/2017MNRAS.467.2377R/abstract).
Further examples of how to work with the SimulationArchive are provided in an [iPython](ipython_examples/SimulationArchive.ipynb) and [C example](c_examples/simulationarchive.md).

## Creating Simulation Archive snapshots

The following code shows how to manually append a Simulation Archive snapshot to a file.
If the file does not exist yet, the function outputs a new binary file. 
If the file already exists, the function will append a Simulation Archive snapshot to the existing file. 

=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... work on simulation ...
    reb_simulationarchive_snapshot(r, "archive.bin");
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.simulationarchive_snapshpt("archive.bin")
    ```
    You can pass the optional argument `deletefile=True` to delete the file if it already exists.
    By default the function appends a snapshot to existing files.

Instead of manually outputting each snapshot, you can also automate this process as shown below.

### Regular time intervals
The following code automatically creates a Simulation Archive snapshot at regular intervals.
== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... work on simulation ...
    reb_simulationarchive_automate_interval(r, "archive.bin", 10.);
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.automateSimulationArchive("archive.bin", interval=10.)
    ```

### Regular number of timesteps
The following code automatically creates a Simulation Archive snapshot after a fixed number of timesteps.
== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... work on simulation ...
    reb_simulationarchive_automate_step(r, "archive.bin", 100);
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.automateSimulationArchive("archive.bin", step=100)
    ```
!!! Info
    This method is in general more reliable than the interval method.
    The reason is that the number of timesteps is an integer valuem whereas the time is a floating point number.
    If you run long simulations, you might encounter issues with finite floating point precision.
    This only affects very high accuracy simulation where you want to make sure the outputs occur exactly at the right timestep. 


### Regular wall-time interbals
The following code automatically creates a Simulation Archive snapshot after a fixed wall-time.
This is particularly useful for creating restart files when running long simulations.
The wall-time is given in seconds.
== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... work on simulation ...
    reb_simulationarchive_automate_walltime(r, "archive.bin", 120.); // 2 minutes between snapshots
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... work on simulation ...
    sim.automateSimulationArchive("archive.bin", walltime=120) # 2 minutes
    ```

## Reading Simulation Archives

### Reading one snapshot
The following example shows how to read in a specific snapshot of a Simulation Archive.
If you pass a negative number for the snapshot, it will wrap around to the end of the Simulation Archive.
For example, the last snapshot in the file would have the index `-1`, the second to last `-2`, and so on.
=== "C"
    ```c
    struct reb_simulationarchive* archive = reb_open_simulationarchive("archive.bin");
    struct reb_simulation* r = reb_create_simulation_from_simulationarchive(archive, 12); // snapshot with index 12
    reb_close_simulationarchive(archive);
    // ... work on simulation ...
    ```

=== "Python"
    ```python
    sim = rebound.Simulation("archive.bin", snapshot=12)
    # ... work on simulation ...
    ```
    In python, `snapshot=-1` is the default. 
    Thus `#!python sim = rebound.Simulation("archive.bin")` will create a new simulation from the last snapshot in the archive. 

