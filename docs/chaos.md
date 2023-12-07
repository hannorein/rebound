# Chaos indicators
REBOUND supports different chaos indicators.
All of these make use of variational equations, but most of the complexity is hidden.

## Initialization
If you want to use a chaos indicator in REBOUND, first add all the particles to your simulations.
Then, initialize the variational particles and MEGNO variables with
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... add particles ...
    reb_simulation_init_megno(r);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    sim.init_megno()
    ```

REBOUND uses random numbers to initialize the variational particles.
The initial seed is chosen based on the current time and the process id. 
This ensures the seed is different every time you run the simulation.
See the discussion on [random sampling](c_randomsamplingfunctions.md) for more details.

If you want to have reproducible result, you can specify the seed manually:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... add particles ...
    reb_simulation_init_megno_seed(r, 0); // 0 is the initial seed
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
    struct reb_simulation* r = reb_simulation_create();
    // ... add particles ...
    reb_simulation_init_megno_seed(r);
    // ... integrate ...
    printf("MEGNO = %f\n", reb_simulation_megno(r));
    printf("LCN = %f\n", reb_simulation_lyapunov(r));
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    sim.init_megno()
    # ... integrate ...
    print("MEGNO", sim.calculate_megno())
    print("LCN", sim.lyapunov())
    ```
!!! Note
    Using chaos indicators is not always straightforward. 
    It can be particularly tricky to figure out how long to integrate for.
    If the integration time is too short, you might not capture the Lyapunov timescale accurately. 
    On the other hand, if the integration time is too long, you can run into problems as well because quantities tend to grow exponentially with time in chaotic systems.

!!! Note
    There are different definitions of the LCN which might differ by a factor of order unity.
    Here, we're following Eq. 24 of [Cincotta and Simo (2000)](https://aas.aanda.org/articles/aas/abs/2000/20/h1686/h1686.html).

## Re-scaling of variational equations 

!!! Important
    This is a new feature, first implemented in version 3.21

REBOUND will automatically re-scale first order variational equation once any coordinate of a variational particle becomes larger than $10^{100}$. 
This is only possible for first order variational equations because they are linear. 
It is not possible to re-scale second order equations. 
For the calculation of MEGNO, all this is done behind the scenes and no user intervention is needed.
However, should you be interested in the actual value of the variational particles, for example to calculate a Lyapunov exponent manually, then you need to take the value of the `lrescale` variable into account. 
This variable contains the natural logarithm of all re-scaling factors that have been applied throughout the integration to a given set of variational particles.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    // ... add particles ...
    reb_simulation_init_megno_seed(r);
    // ... integrate ...
    struct reb_variational_configuration* vc = &(r->var_config[0]);
    double log_x = log(r->particles[vc->index].x) + vc->lrescale; // log of x coordinate of variational particle
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    # ... add particles ...
    sim.init_megno()
    # ... integrate ...
    vc = sim.var_config[0]
    log_x = vc.particles[0] + vc.lrescale // log of x coordinate of variational particle

Note that above the  calculations involving the rescaling factors have been done in log space as floating point numbers cannot be used to represent a number larger than $10^{308}$.
For a more complete example, check out the [iPython Variational Equation example](ipython_examples/VariationalEquations.ipynb).
