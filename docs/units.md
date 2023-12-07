# Units
## Using G=1
By default, REBOUND simulations use units in which $G=1$.
That might be confusing at first if you're used to working in SI units.
The reason for setting $G=1$ is that gravity is scale-free.
Imagine a simulation with two particles orbiting each other.
For REBOUND, it doesn't matter if this is a planet orbiting a star, a moon orbiting a planet, or a spacecraft orbiting a moon.
It is only a matter of interpreting the simulation. 

As an example, suppose we use $G=1$, have a central object of mass $M$, and a test particle orbiting on a circular orbit at a distance $a=1$. 
This scenario can be setup with the following code:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt("m", 1.);
    reb_simulation_add_fmt("a", 1.);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1.)
    sim.add(a=1.)
    ```
The orbital period of this binary is given by 

$$
P = 2\pi\sqrt{\frac{a^3}{GM}} = 2\pi
$$

We can confirm this by calculating the orbital period with REBOUND:
=== "C"
    ```c
    struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);
    printf("P=%f\n", o.P);
    ```
=== "Python"
    ```python
    print(sim.particles[1].P)
    ```
If we interpret the central object as the sun, and the test particle as the Earth, then $M=1$ corresponds to one solar mass and $a=1$ corresponds to one astronomical unit. 
We know that the Earth takes one year for one orbit around the sun. 
Thus, one year corresponds to $2\pi$ in these units.

An alternative interpretation of the same REBOUND simulation could be the following. 
Suppose the central object is the Earth. 
Then $M=1$ corresponds to one Earth mass.
If we consider the test particle to be the International Space Station, then $a=1$ corresponds to $6790{\rm km}$ ($420{\rm km}$ above MSL).
The orbit of the particle still has a period of $2\pi$ in our units, but this would now correspond to 92.8 minutes.

!!! Info
    One advantage of keeping $G=1$ is that you can choose units where all number in REBOUND have a magnitude of around one, as in the example above. 
    If you where to choose other units involving centimetres or seconds, then you would have to deal with very large or very small numbers.

## Changing G
If you prefer to change the value of $G$, you can!
The following example sets $G$ to its value in SI units:
=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    r->G = 6.6743e-11; //  m^3 / kg s^2
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.G = 6.6743e-11 # m^3 / kg s^2
    ```
From now on, all quantities that have unit of length (semi-major axis, particle radius, etc) need to be specified (and will be output) in meters.
All quantities that have units of time (timestep, orbital period, etc) need to be specified in seconds. 
All quantities that have units of mass need to be specified in kg.
All quantities that have units of velocity need to be specified in meters per second.
And so on.

## Convenience methods in python
The python version of REBOUND comes with its own set of convenience functions for changing the system of units.
Check out the [iPython example](ipython_examples/Units.ipynb).
