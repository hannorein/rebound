# Chaos indicators
REBOUND supports different chaos indicators.
All of these make use of variational equations, but most of the complexity is hidden.

## Initialization
If you want to use a chaos indicator in REBOUND, first add all the particles to your simulations.
Then, initialize the variational particles and MEGNO variables with
=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... add particles ...
    reb_tools_megno_init(r);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    sim.init_megno()
    ```

REBOUND uses random numbers to initialize the variational particles.
The initial seed is chosen based on the current time and the process id. 
This ensures the seed if different every time you run the simulation.
See the discussion on [random sampling](c_randomsamplingfunctions.md) for more details.

If you want to have reproducible result, you can specify the seed manually:
=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... add particles ...
    reb_tools_megno_init_seed(r, 0); // 0 is the initial seed
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    sim.init_megno(seed=0) # 0 is the initial seed
    ```
## Accessing chaos indicators
Once you've initialized the chaos indicators, you can integrate the simulation normally.
To print out the MEGNO value or the largest Lyapunov characteristic number (LCN), use the following syntax:

=== "C"
    ```c
    struct reb_simulation* r = reb_create_simulation();
    // ... add particles ...
    reb_tools_megno_init_seed(r);
    // ... integrate ...
    printf("MEGNO = %f\n", reb_tools_calculate_megno(r));
    printf("LCN = %f\n", reb_tools_calculate_lyapunov(r));
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    sim.init_megno()
    # ... integrate ...
    print("MEGNO", sim.calculate_megno())
    print("LCN", sim.calculate_lyapunov())
    ```
!!! Important
    Using chaos indicators is not always straightforward. 
    It can be particularly tricky to figure out how long to integrate for.
    If the integration time is too short, you might not capture the Lyapunov timescale accurately. 
    On the other hand, if the integration time is too long, you can run into problems as well because quantities tend to grow exponentially with time in chaotic systems.
    In the worst case, the magnitude of the variational particles' coordinates will reach the maximum range of double precision floating point numbers. 
    Ideally, you would want to rescale the variational particles whenever that happens. 
    However, this is currently not implemented in these convenience methods and you will have to use variational particles directly.
