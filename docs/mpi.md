# MPI (Message Passing Interface)

!!! info inline end Python
    Only the C version of REBOUND supports MPI. The python version of REBOUND does not support MPI.


REBOUND supports parallelization on distributed memory systems with MPI (Message Passing Interface). 
This can be useful to accelerate simulations but it only makes sense if certain conditions are met:

- A large number of particles (at least a few thousand) is needed for the parallelization to provide a speed-up. If the number of particles is too small, then the parallelization will slow down the simulation because the communication will be the new bottleneck. 
- Only simulations that use a tree code can be parallelized. The tree is used for the domain decomposition.


Use cases where MPI might be a good way to speed up simulations are:

- [A self-gravitating disk](../c_examples/selfgravity_disc_mpi/)
- [A shearing sheet simulation](../c_examples/shearing_sheet_mpi/) of self-gravitating or collisional particles (e.g. to simulate Saturn's Rings)

## Basic Workflow
The basic workflow when using MPI is as follows. You need to enable MPI and choose the appropriate compiler for MPI in the Makefile of the problem directory:

```
export MPI=1 
export CC=mpicc
```

After creating the simulation structure, you need to initialize the tree structure:

``` c
struct reb_simulation* r = reb_simulation_create();
r->gravity          = REB_GRAVITY_TREE;
r->collision        = REB_COLLISION_TREE;
// other configuration options
reb_simulation_configure_box(r, boxsize, 2, 2, 1);
```    

The number of root trees needs to be an integer multiple of the number of MPI processes.
In this example, 2 root boxes are used in the x and y directions and 1 in the z direction.
Thus you can use 1, 2, or 4 MPI processes.

The combined size of the trees should be large enough to contain all your particles.
After the tree has been initialized, you need to initialize MPI:

``` c
reb_mpi_init(r);
```

You can now add particles and integrate the simulation. 
Once the simulation is done terminates the MPI execution environment and cleanup the memory used by the simulation with:

``` c
reb_mpi_finalize(r);
reb_simulation_free(r);
```

How to submit and run parallel jobs depends on your computing cluster. Please contact your cluster administrator if you have questions about this.

## Support
In general, using REBOUND with MPI requires a lot more work on the user's side to make thing work. 
Many features are currently not compatible with MPI, or require some extra thought, for example binary input/output and Simulationarchives.
If you would like to use a features with MPI that is currently not supported, or you have any other questions regarding MPI and REBOUND, please [open an issue on GitHub](https://github.com/hannorein/rebound/issues).

