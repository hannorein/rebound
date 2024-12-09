# Operators

## Copying simulations
REBOUND makes it very easy to copy a simulation. 
This can be very helpful in many cases.
For example, you can keep a record of your initial conditions by simply making a copy of the simulation before you start any integration. 


=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    struct reb_simulation* r_copy = reb_simulation_copy(r);
    ```
=== "Python"
    ```python
    r = rebound.Simulation()
    r_copy = r.copy()
    ```
!!! Info
    The above function calls create a deep copy of the simulation.
    All the data in the simulation is duplicated, including the particle data.
    If you use function pointer in the original simulation, you will need to manually reset them.

## Adding, subtracting, multiplying simulations
REBOUND allows you to manipulate entire simulations with 'arithmetic' operations.
For example:

=== "C"
    ```c
    struct reb_simulation* r1 = reb_simulation_create();
    struct reb_simulation* r2 = reb_simulation_create();
    // ... setup simulations ...
    reb_simulation_isub(r1, r2);
    reb_simulation_iadd(r1, r2);
    ```
=== "Python"
    ```python
    r1 = rebound.Simulation()
    r2 = rebound.Simulation()
    # ... setup simulations ...
    r1 -= r2
    r1 += r2
    ```

In the above example, each particle in `r2` is subtracted from the corresponding particle in `r1`, in the sense described above (element-wise position, velocity and mass), then the operation is reversed in the next line when the simulation are added together.  
This feature can come in very handy when you want to compare two different simulations. 
For example, you can run two simulations with different timesteps and then subtract the simulation from each other after the integration to how much of a difference the timestep makes.
These operations will fail if the number of particles are not the same in `r1` and `r2`.

You can also multiply a simulation with two scalars, one for the position and one for the velocity coordinates.
This can come in handy when re-scaling a simulation

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation ...
    reb_simulation_imul(r, 2., 3.,);
    ```
=== "Python"
    ```python
    r = rebound.Simulation()
    # ... setup simulation ...
    r.multiply(2., 3.)
    ```
In the above the position coordinates of all particles are multiplied by 2, all velocity coordinates are multiplied by 3.

## Comparing simulations
You can compare if simulations are equal to each other using the following syntax:
=== "C"
    ```c
    struct reb_simulation* r1 = reb_simulation_create();
    struct reb_simulation* r2 = reb_simulation_create();
    // ... setup simulations ...
    if (reb_simulation_diff(r1, r2, 2)){
        // Simulations are NOT equal
    }
    ```
    For debugging purposes, it can be useful to print out the differences. 
    This is done by passing `1` as the last argument:
    ```c
    struct reb_simulation* r1 = reb_simulation_create();
    struct reb_simulation* r2 = reb_simulation_create();
    // ... setup simulations ...
    reb_simulation_diff(r1, r2, 1); // prints out diferences
    ```
=== "Python"
    ```python
    r1 = rebound.Simulation()
    r2 = rebound.Simulation()
    # ... setup simulations ...
    if r1 == r2:
        print("Simulations are equal")
    ```
