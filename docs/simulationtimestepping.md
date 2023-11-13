# Timestepping

## Integrate 
In most cases, you will want to integrate your simulation to a given time. 
This could be the time at which you want to create the next output, or a very long time into the future, if you are waiting for an exception to happen (close encounter, ejection, etc). 

In those cases use this syntax:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation, set timestep ...
    reb_simulation_integrate(r, 100.); // integrate until t=100.
    ```
    If you want to integrate indefinitely, you can use
    ```c
    reb_simulation_integrate(r, INFINITY); 
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.integrate(100.) # integrate until t=100.
    ```

The integrate function will integrate the simulation until it reaches exactly the time requested.
In most cases the time requested will not be an exact multiple of the timestep, so the timestep will have to be reduced during the last timestep. 
After the requested time has been reached, the timestep will be reverted back to its original value. 

There are cases where you don't want to reduce the timestep, for example in long term integrations with symplectic integrators.
In those cases, you can ask REBOUND to integrate up to a given time and overshoot the requested time by a fraction of the timestep.
This allows REBOUND to maintain a constant timestep throughout the integration.
The following code shows you how to do that.
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation, set timestep ...
    r->exact_finish_time = 0;
    reb_simulation_integrate(r, 100.); // integrate until t=100. or a bit further
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... setup simulation, set timestep ...
    sim.integrate(100., exact_finish_time=0) # integrate until t=100. or a bit further
    ```

If you want to stop a current integration after the current timestep, for example from within the heartbeat function, you can call:
=== "C"
    ```c
    reb_simulation_stop(r); 
    ```
=== "Python"
    sim.stop()
    ```

Note that you might need to manually synchronize the simulation afterwards if you have safe mode turned off.


## Single step

Rather than integrating up to a fixed time, you can also advance the simulation by a single timestep:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation, set timestep ...
    reb_simulation_step(r);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... setup simulation, set timestep ...
    sim.step()
    ```

## Multiple steps
And finally, you can ask REBOUND to advance the simulation by a finite number of steps.
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation, set timestep ...
    reb_simulation_steps(r, 100); // 100 steps
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... setup simulation, set timestep ...
    sim.steps(100) # 100 steps
    ```

## Synchronizing
Depending on the `safe_mode` flag, some integrators perform optimizations which effectively leave a timestep unfinished.
You can manually 'synchronize' the simulation by calling

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation ...
    reb_synchronize(r); 
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... setup simulation ...
    sim.synchronize()
    ```

See the [discussion on integrators](integrators.md) for more information about the `safe_mode` and synchronizing simulations.
