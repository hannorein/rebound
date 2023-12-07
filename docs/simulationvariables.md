# Variables
The following example shows how to access variables in the simulation structure.  
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->G = 1.0;             // Set the gravitational constant
    printf("%f\n", r->t);   // print current simulation time
    ```

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.G = 1.0    # Set the gravitational constant
    print(sim.t)   # print current simulation time
    ```

Below, we list the important variables in the simulation structure. 
To keep the documentation concise, variables which are only intended for internal use are not documented here. 

## Timestepping

`#!c double t`                  
:   Current simulation time. The default value is 0. The value increases if a simulation is integrated forward in time ($dt>0$). See also the [discussion on units](units.md).

`#!c double dt`                 
:   This is the current timestep. The default is 0.001. 
    Make sure to set the timestep to a small fraction (a few percent) of the shortest dynamical timescale in the problem.
    Adaptive integrators such as [IAS15](../integrators/#ias15) will use this value as their initial guess during the first timestep.
    In subsequent timesteps, adaptive integrators will change this value.
    See also the [discussion on units](units.md).
     
`#!c double dt_last_done`       
:   REBOUND sets this variable to the last timestep used. Do not set this variable manually.

`#!c unsigned long long steps_done` 
:   Number of timesteps completed.

`#!c int exact_finish_time`     
:   If this variable is set to 1 (default), then REBOUND will integrate the simulation exactly up to the requested time. 
    Unless the requested time is a multiple of the timestep, REBOUND will need to reduce the timestep to achieve this.
    Set this variable to 0 and REBOUND will not reduce the timestep and will instead overshoot the integration by a fraction of one timestep. 

    === "C"
        ```c
        struct reb_simulation* r = reb_simulation_create();
        r->integrator = REB_INTEGRATOR_LEAPFROG; // uses fixed timestep
        r->dt = 10;
        r->exact_finish_time = 0;
        reb_simulation_integrate(r, 15);
        printf("%f\n", r->t);     // will print 20
        
        r->exact_finish_time = 1; // default
        reb_simulation_integrate(r, 25);
        printf("%f\n", r->t);     // will print 25
        ```

    === "Python"
        In python, you do not need to set this flag in the simulation structure. 
        Instead, you pass it as an argument when calling `integrate()`:
        ```python
        sim = rebound.Simulation()
        sim.integrator = "leapfrog" # uses fixed timestep
        sim.dt = 10
        sim.integrate(15, exact_finish_time=0)
        print(sim.t)                # will print 20

        sim.integrate(25, exact_finish_time=1)
        print(sim.t)                # will print 25
        ```

`#!c double walltime`           
:   This variable keeps track of the wall-time (in seconds) used by REBOUND for this simulation.
    This is counting only the integration itself and not the visualization, heartbeat function, etc.

`#!c void (*heartbeat) (struct reb_simulation* r)`
:   The `heartbeat` function pointer is called at the beginning of the simulation and at the end of each timestep.
    You can use this function to keep track of your simulation, terminate it, or output data.

    === "C"
        ```c
        void heartbeat(struct reb_simulation* r){
            printf("%f\n",r->t);
        }
        int main(int argc, char* argv[]) {
            struct reb_simulation* r = reb_simulation_create();
            r->heartbeat = heartbeat;
            // ...
        }
        ```
    
    === "Python"
        ```python
        def heartbeat(sim_pointer):
            sim = sim_pointer.contents
            print(sim.t)
        
        sim = rebound.Simulation()
        sim.heartbeat = heartbeat
        # ...
        ```

`#!c void (*pre_timestep_modifications) (struct reb_simulation* const r)`

`#!c void (*post_timestep_modifications) (struct reb_simulation* const r)`
:   Similar to the heartbeat function, these function pointers allow you to make changes before and after each timestep.
    These pointers are also used by REBOUNDx.

## Gravity

`#!c double G`                  
:   Gravitational constant. By default, this value is 1. 
    If $G=1$, then an orbit with semi-major axis $a=1$ has a period of $P=2\pi$.
    See also the [discussion on units](units.md).

`#!c double softening`          
:   This is the gravitational softening parameter. 
    The gravitational force of a particle in the $x$ direction is calculated as
    $F_x = -x \frac{G m_1 m_2}{(x^2 +y^2 +z^2 + b^2)^{3/2}}$, where $b$ is the gravitational softening parameter.
    This can be used to remove strong force gradients on small scales, e.g. during close encounters.
    The default is 0 (no softening). 

`#!c double opening_angle2`     
:   This variable determines the accuracy of the gravity calculation when the tree bases gravity routine is used.
    It is the square of the cell opening angle $\theta$. 
    See [Rein & Liu](https://ui.adsabs.harvard.edu/abs/2012A%26A...537A.128R/abstract) for a discussion of the tree code.

`#!c unsigned int force_is_velocity_dependent` 
:   If this variable is set to 0 (default), then the force can not contain velocity dependent terms.
    Setting this to 1 is slower but allows for velocity dependent forces (e.g. drag force). 
    Note that gravitational forces alone are not velocity dependent. 
    
`#!c unsigned int gravity_ignore_terms`        
:   This variable determines if the gravity form the central object is included in the gravity calculation.
    In general the integrators will set this variable automatically and nothing needs to be changed by the user.
    Possible values are:

    - 0 include all terms
    - 1 ignore terms not required for WHFast with Jacobi coordinates
    - 2 ignore terms not required for WHFast with democratic heliocentirc coordinates

`#!c void (*additional_forces) (struct reb_simulation* const r)`
:   This function allows the user to add additional (non-gravitational) forces.
     
    !!! Todo
        Add examples.

## Particles

`#!c struct reb_particle* particles` 
:   All particles are stored in this array.
    A particle is represented by the `reb_particle` structure in C.
    The python class `Particle` is an abstraction of the `reb_particle` structure in C.
    
    The following example changes a particle's x coordinate:
    
    === "C"
        ```c
        struct reb_simulation* r = reb_simulation_create();
        // ... setup simulation, add particles ...
        r->particles[0].x = 1;
        ```

    === "Python"
        ```python
        sim = rebound.Simulation()
        # ... setup simulation, add particles ... 
        sim.particles[0].x = 1
        ```
    
    The memory for this array is managed by REBOUND.
    To add and remove particles, don't modify this array directly. 
    Instead use the `_add` and `_remove` functions.

    The order in which particles are added matters in multiple situation:

    - When test particles are used, active particles need to be added before test particles.
    - When integrating systems with the WHFast integrator, the central object needs to be added first.
    - When Jacobi coordinates are used, then the particles needs to be added from the inside out (star, inner planet, outer planet).

`#!c int N`                     
:   Current number of particles in this REBOUND simulation. 
    This number includes all active, test, and variational particles.
    The simulation stops when this number is 0 and there are no more particles in the simulation.
    The default is 0.

`#!c int N_active`              
:   This is the number of active particles in the simulation.
    Only active particles contribute to the force in the gravity calculation.
    The default is -1 which means the number of active particles is equal to the number of particles, `N`. 
    Particles with an index larger or equal than `N_active` are considered test-particles.
   
    The following example sets up a simulation with two active particles and one massless test-particle.
    === "C"
        ```c
        struct reb_simulation* r = reb_simulation_create();
        reb_simulation_add_fmt(r, "m", 1.0);
        reb_simulation_add_fmt(r, "m a", 1e-3, 1.0);
        reb_simulation_add_fmt(r, "m a", 0.0, 2.0);
        r->N_active = 2;
        ```

    === "Python"
        ```python
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3, a=1)
        sim.add(m=0, a=2)
        sim.N_active = 2
        ```

`#!c int testparticle_type`     
:   This determines the type of the particles with `index >= N_active`. 
    REBOUND supports two different test-particle types:

    - If this variable is set to 0, then test particle does not influence any other particle (default).
    - If this variable is set to 1, then active particles (those with `index < N_active`) feel test-particles (similar to MERCURY's small particles). 

    Test-particles never feel each other.

`#!c int N_var`                 
:   Total number of variational particles. Default: 0.

`#!c int N_var_config`          
:   Number of variational particle configurations. Default: 0.


## Collisions


`#!c int (*collision_resolve) (struct reb_simulation* const r, struct reb_collision)` 
:   This is a function pointer which determines how a collision is resolved. By default, it is NULL, assuming hard sphere model.
    A return value of 0 indicates that both particles remain in the simulation. A return value of 1 (2) indicates that particle 1 (2) should be removed from the simulation. A return value of 3 indicates that both particles should be removed from the simulation. 
    See [the discussion on collisions](collisions.md#resolving-collisions) for more information on how to use this function pointer. 

`#!c int track_energy_offset`   
:   Set this variable to 1 to track energy change during collisions and ejections (default: 0).
    This is helpful if you want to keep track of an integrator's accuracy and physical collisions do not conserve energy.

`#!c double energy_offset`      
:   Energy offset due to collisions and ejections (only calculated if `track_energy_offset=1`).

`#!c int collision_resolve_keep_sorted` 
:   If set to 1, then particles are kept sorted when a particle is removed during a collision.

`#!c double minimum_collision_velocity`  
:   When collisions are resolved with the hard sphere collision resolve function, then the post impact velocity between the two particles will be at least as large as this value. Default 0. Setting this to a value larger than zero might prevent particles sinking into each other. 

`#!c double collisions_plog`    
:   This variable keeps track of momentum exchange during collisions. This can be used to calculate collisional viscosity in ring systems.

`#!c long collisions_log_n`      
:   This variable keeps track of the number of collisions that have occurred. This can be used to calculate statistical quantities of collisional systems.

`#!c double (*coefficient_of_restitution) (const struct reb_simulation* const r, double v)`
:   This is a callback function which gets called when a hard-sphere collision occurs and the coefficient of restitution is required.
    By default, this function pointer is NULL and a coefficient of restitution of 1 is assumed.
    The impact velocity of the collision is given to allow for velocity dependent coefficients of restitution.

## Miscellaneous 

`#!c enum REB_STATUS status`    
:   This variable indicates the current status of the simulation. By setting this to 1, one can force a graceful exit at the end of the next timestep.    

`#!c double exit_max_distance`  
:   The integration will stop if any particle is further away from origin than this value.

`#!c double exit_min_distance`  
:   The integration will stop if any two particles come closer together than this value.

`#!c double usleep`             
:   Sleep this number of microseconds after each timestep. This can be useful for slowing down the simulation, for example for rendering visualizations.  

`#!c int N_ghost_x, N_ghost_y,  N_ghost_z`               
:   Number of ghost-boxes in x, y, and z directions. 

`#!c unsigned int rand_seed`    
:   Seed for random number generators. This will be automatically initialized to a random number based on the current time and the process id. However, it can also be set manually to make the simulation reproducible and always return the same sequence of random numbers.


## Module selection 

The following variables in the simulation structure determine which modules are selected. 
The [gravity solvers](gravity.md), [collision detection algorithms](collisions.md), [boundary conditions](boundaryconditions.md), and [integration methods](integrators.md) are explained in detail on their own pages.

`#!c enum visualization`

`#!c enum collision`

`#!c enum integrator`

`#!c enum boundary`

`#!c enum gravity`

## Integrator configuration 

The following variables in the simulation structure contain the configuration for the individual integrators. 
They are described on their own [separate page](integrators.md). 

`#!c struct reb_integrator_sei ri_sei`

`#!c struct reb_integrator_whfast ri_whfast`

`#!c struct reb_integrator_saba ri_saba`

`#!c struct reb_integrator_ias15 ri_ias15`

`#!c struct reb_integrator_mercurius ri_mercurius`

`#!c struct reb_integrator_janus ri_janus`

`#!c struct reb_integrator_eos ri_eos`


