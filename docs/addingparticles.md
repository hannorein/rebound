# Adding particles

![type:video](https://www.youtube.com/embed/FoTwDtAeJyk)

Once you've created a [simulation object](simulation.md), you can add particles to it. 
REBOUND supports several different ways to do that. 
Also check out the [discussion on particle operators](particleoperators.md).

## Adding particles manually
One way to add a particle to a simulation is to first manually create a particle object, then calling a function to add the particle to the simulation.
Because the function will make a copy of the particle, you can safely delete the original particle object after you've added it to a simulation. 
The following code shows an example on how to add particles this way:

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle p = {0};
    p.m = 1.;
    p.x = 1.;
    reb_simulation_add(r, p);
    ```
    !!! Important
        The `= {0}` syntax above ensures that the struct is initialized with zeros.
        Otherwise, you need to set every member of the struct to ensure that there are no
        uninitialized values.


=== "Python"
    ```python
    sim = rebound.Simulation()
    p = rebound.Particle()
    p.m = 1.
    p.x = 1.
    sim.add(p)
    ```

You can also use orbital parameters to initialize the particle object.
=== "C"
    In C, this is done by calling the `reb_particle_from_orbit` function. Its arguments are gravitational constant, primary object, mass, semi-major axis, eccentricity, inclination, longitude of ascending node, argument of pericenter, and true anomaly. 
    It returns an initialized particle object which you can then add to the simulation.
    ```c
    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle primary = {0};
    primary.m = 1;
    reb_simulation_add(r, primary);
    struct reb_particle planet = reb_particle_from_orbit(r->G, primary, 1e-3, 1., 0., 0., 0., 0., 0.);
    reb_simulation_add(r, planet);
    ```

    You can also the coordinates described by [Pal 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.396.1737P/abstract) to initialize orbits using the following function:

    ```c
    struct reb_particle reb_particle_from_pal(double G, struct reb_particle primary, double m, double a, double lambda, double k, double h, double ix, double iy);
    ``` 
    Here, `lambda` is the longitude, `h` is $e\cos(\omega)$, `k` is $e\sin(\omega)$, `ix` and `iy` are the x and y components of the inclination respectively. 

=== "Python"
    In python, you can create and initialize particles using the constructor of the `Particle` class.
    ```python
    sim = rebound.Simulation()
    primary = rebound.Particle(m=1., x=1.)
    sim.add(primary)
    ```
    If you want to use orbital parameters, you need to pass the primary and the simulation to the constructor:
    ```python
    planet = rebound.Particle(simulation=sim, primary=primary, m=1e-3, a=1., e=0.1)
    ```
    You can use any combination of orbital parameters that makes physically sense. 
    See [the discussion on orbital elements](orbitalelements.md) for more details.

    !!! Note
        In most cases you can simply use the convience function described below.
        This way you don't have to create a particle object just to add it to the simulation.

## Convenience functions
By far the easiest way to add particles to REBOUND is to use a convenience function.
=== "C"
    In C, the function is called `reb_simulation_add_fmt` and has the following syntax:
    ```c
    void reb_simulation_add_fmt(struct reb_simulation* r, const char* fmt, ...);
    ```
    This is a [variadic function](https://en.cppreference.com/w/c/variadic) which takes a variable number of arguments similar to the `printf` function.
    The following code shows how this function is used.
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m", 1.0);                // star at origin with mass 1
    reb_simulation_add_fmt(r, "m a", 1e-3, 1.0);        // planet with mass 1e-3 and semi-major axis 1
    reb_simulation_add_fmt(r, "m a e", 1e-3, 2.0, 0.1); // planet with mass 1e-3, semi-major axis 2, and eccentricity 0.1
    reb_simulation_add_fmt(r, "m x vy", 1e-6, 1., 1.);  // planet with mass 1e-6, cartesian coordinates
    ```

    The first argument is the simulation to which you want to add the particle.
    The second argument is a format string and it determines how many other arguments the function expects.

    !!! Danger
        You need to pass exactly the right number of arguments to `reb_simulation_add_fmt` as indicated by your format string. 
        Each argument also has to be the right type (mostly double floating point numbers).
        The latter is particularly important. If you call the function like this:
        ```c
        reb_simulation_add_fmt(r, "m a", 1, 1);
        ```
        then the arguments are integers, not doubles. This can lead to unexpected behaviour that is very difficult to debug.
        The correct way to call the function is by making sure the arguments are doubles (by adding a `.`):
        ```c
        reb_simulation_add_fmt(r, "m a", 1.0, 1.0);
        ```

    The following parameters are supported:

    Parameter | Description
    --------- | -----------
    `m`|          mass  (default: 0)
    `x, y, z`|    positions in Cartesian coordinates  (default: 0)
    `vx, vy, vz`| velocities in Cartesian coordinates (default: 0)
    `primary`|    primary body for converting orbital elements to cartesian (default: center of mass of the particles in the passed simulation, i.e., this will yield Jacobi coordinates as one progressively adds particles) 
    `a`|          semi-major axis (a or P required if passing orbital elements)
    `P`|          orbital period (a or P required if passing orbital elements)
    `e`|          eccentricity                (default: 0)
    `inc`|        inclination                 (default: 0)
    `Omega`|      longitude of ascending node (default: 0)
    `omega`|      argument of pericenter      (default: 0)
    `pomega`|     longitude of pericenter     (default: 0)
    `f`|          true anomaly                (default: 0)
    `M`|          mean anomaly                (default: 0)
    `E`|          eccentric anomaly           (default: 0)
    `l`|          mean longitude              (default: 0)
    `theta`|      true longitude              (default: 0)
    `T`|          time of pericenter passage  
    `h, k, ix, iy`|      See [Pal 2009](https://ui.adsabs.harvard.edu/abs/2009MNRAS.396.1737P/abstract) for a definition  (default: 0)
    `r`|          physical particle radius
 
    You can use any combination of these parameters at the same time. 
    If a combination is unphysical, no particle will be added and an error will be outputted.
    For example, you can only specify one longitude or anomaly.

=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1)                 # star at origin with mass 1                                      
    sim.add(m=1e-3, a=1.)        # planet with mass 1e-3 and semi-major axis 1
    sim.add(m=1e-3, a=2., e=0.1) # planet with mass 1e-3, semi-major axis 2, and eccentricity 0.1
    sim.add(m=1e-6, x=1., vy=1.) # planet with mass 1e-6, cartesian coordinates
    ```

See [the discussion on orbital elements](orbitalelements.md) for more details.


## Solar System planets
If you want to quickly try something out, you can use a set of initial conditions for the Solar System that come with REBOUND: 

```python
sim = rebound.Simulation()
rebound.data.add_solar_system(sim)
```

and similarly for the outer Solar System:

```python
sim = rebound.Simulation()
rebound.data.add_outer_solar_system(sim)
```


This is currently only supported in python.

!!! Note
    These initial conditions are intended for testing integration methods. They might not be very accurate and should not be used for detailed dynamical studies of the Solar System.


