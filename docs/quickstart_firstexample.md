# Your first REBOUND simulation

If you have successfully [installed REBOUND](quickstart_installation.md), then you are now ready to run your first simulation.
On this page, we'll walk you through one simple example, line by line.

## Python version

There are different ways to run python code:

- interactively, by using the python interpreter
- by executing a python script
- by using a Jupyter notebook

All of these methods work with REBOUND.
Choose whichever you are most comfortable with.

We start by importing REBOUND:

```python
import rebound
```

To run an N-body simulation, we need to create a simulation object first:

```python
sim = rebound.Simulation()
```

Then, we [add particles](addingparticles.md) to the simulation:

```python
sim.add(m=1.)                # Central object
sim.add(m=1e-3, a=1., e=0.1) # Jupiter mass planet 
sim.add(a=1.4, e=0.1)        # Massless test particle
```

We are working in units where $G=1$. [Click here](units.md#using-g1) to learn more about what these units mean. 
Now we can integrate the particles forward in time using the default integrator ([IAS15](integrators.md#ias15)) for 100 time units:

```python
sim.integrate(100.)
```

Finally, let us output the Cartesian coordinates and the orbital parameters at the end of the simulation:

```python
for p in sim.particles:
    print(p.x, p.y, p.z)
for o in sim.orbits(): 
    print(o)
```

As a next step, have a look at the examples and tutorials in the `python_examples` and `ipython_examples` directories.

## C version

A very short example is provided in the `examples/simplest/` directory. 
Go to this directory with 

```bash
cd examples/simplest/
```

Then have a look at the source code in the `problem.c` file. First, we include the REBOUND header file which contains all the public function prototype and datatype definitions for REBOUND:

```c
#include "rebound.h"
```

In the main function, we first create a REBOUND simulation with 

```c
struct reb_simulation* r = reb_simulation_create();
```

This function has now allocated memory for the simulation and initialized all the variables in the simulation to their default values.
We can then [add particles](addingparticles.md) to the simulation:

```c
reb_simulation_add_fmt(r, "m", 1.);                // Central object
reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
reb_simulation_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle
```

We are working in units where $G=1$. [Click here](units.md#using-g1) to learn more about what these units mean. 
We then integrate the simulation for 100 time units with the default integrator ([IAS15](integrators.md#ias15)):

```c
reb_simulation_integrate(r,100.);
```

After the integration is done, we can output the Cartesian coordinates and the orbital parameters:

```c
for (int i=0; i<r->N; i++){
    struct reb_particle p = r->particles[i];
    printf("%f %f %f\n", p.x, p.y, p.z);
}
struct reb_particle primary = r->particles[0];
for (int i=1; i<r->N; i++){
    struct reb_particle p = r->particles[i];
    struct reb_orbit o = reb_orbit_from_particle(r->G, p, primary);
    printf("%f %f %f\n", o.a, o.e, o.f);
}
```

To compile the example, simple type

```bash
make
```

into a terminal window while you're in the `examples/simplest/` directory. Then run the simulation with 

=== "Linux/Mac"
    ```bash
    ./rebound

=== "Windows"
    ```bash
    rebound.exe
    ```
