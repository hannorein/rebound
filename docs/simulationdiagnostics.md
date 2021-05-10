# Diagnostics 

## Energy
You can calculate the total energy (kinetic plus potential energy) of a simulation using the following function:

=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    double energy = reb_calculate_energy(r);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    energy = sim.calculate_energy()
    ```

## Angular momentum
You can calculate the angular momentum of a simulation using the following function:

=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    double energy = reb_calculate_angular_momentum(r);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    Lx, Ly, Lz = sim.calculate_angular_momentum()
    ```

## Center-of-mass
You can calculate the center-of-mass of a simulation using the functions below. 
The return value is particle object with mass, position, and velocity reflecting those of the center-of-mass.

=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... setup simulation ...
    struct reb_particle com = reb_get_com(r);
    ```
    You can also return the center-of-mass for particles with indices in a given range.
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... setup simulation ...
    struct reb_particle com = reb_get_com_range(r, 6, 9);
    ```
    In the above example, the particles 6, 7, and 8 are included in the calculation.

=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... setup simulation ...
    com = sim.calculate_com()
    ```

