# Moving reference frames
Compared to other N-body codes, REBOUND does not use a predefined coordinate system. 
It works in any inertial frame.
This makes setting up simulations and interpreting outputs more intuitive. 
However, one often wants to move to a specific coordinate system.
REBOUND has several built-in functions to do that.

## Heliocentric frame
Here is how you can move the simulation to the heliocentric frame.
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation ...
    reb_simulation_move_to_hel(r);
    ```
=== "Python"
    ```python
    r = rebound.Simulation()
    # ... setup simulation ...
    r.move_to_hel()
    ```
This moves all particles in the simulation by the same amount so that afterwards, the particle with index 0 is located at the origin. 
Note that as the integration progresses, it is not guaranteed that the particle with index 0 remains at the origin. 
Most likely, it will drift away from the origin. 
Therefore, if you require outputs in the heliocentric frame, call `move_to_hel` before you create an output. 
Variational equations are not affected by this operation.

## Center-of-mass frame
You can also move a simulation to the center-of-mass frame, the inertial frame where the center-of-mass is at the origin. 
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... setup simulation ...
    reb_simulation_move_to_com(r);
    ```
=== "Python"
    ```python
    r = rebound.Simulation()
    # ... setup simulation ...
    r.move_to_com()
    ```
!!! Important
    If you are not in the center-of-mass frame, the center-of-mass will and all the particles will slowly drift away from the origin. 
    This has important consequences for long-term integrations. 
    If the particles are far away from the origin, you might increase the numerical errors due to finite floating point precision. 
    By moving to the center-of-mass frame after setting up all the particles, you avoid these issues.

