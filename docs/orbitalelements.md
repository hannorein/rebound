# Orbital elements

This section discusses orbital parameters.
We focus on the implementation and conventions in REBOUND.


The following image illustrated the most important angles used.
In REBOUND the reference direction is the positive x direction, the reference plane
is the xy plane.

![Orbital elements](img/orbit.png "Image from wikipedia. CC-BY-SA-3.")

## Orbit structure 

Variable name   | Description
--------------- | ------------
`d`             | radial distance from reference 
`v`             | velocity relative to central object's velocity
`h`             | specific angular momentum
`P`             | orbital period (negative if hyperbolic)
`n`             | mean motion    (negative if hyperbolic)
`a`             | semimajor axis
`e`             | eccentricity
`inc`           | inclination
`Omega`         | longitude of ascending node
`omega`         | argument of pericenter
`pomega`        | longitude of pericenter
`f`             | true anomaly
`M`             | mean anomaly
`E`             | Eccentric anomaly. Because this requires solving Kepler's equation it is only calculated when needed in python and never calculated in C. To get the eccentric anomaly in C, use the function `double reb_tools_M_to_E(double e, double M)`
`l`             | mean longitude = Omega + omega + M
`theta`         | true longitude = Omega + omega + f
`T`             | time of pericenter passage
`rhill`         | Hill radius, $r_{\rm hill} =a\sqrt[3]{\frac{m}{3M}}$

!!! Important
    All angles in REBOUND are in radians. 
    Variables which have length, time or velocity units use code units.

### Particle to orbit

The following function allows you to calculate the orbital elements of a particle.

=== "C"
    ```c
    struct reb_simulation* r = create_simulation();
    reb_add_fmt(r, "m", 1.); // star
    reb_add_fmt(r, "a e", 1., 0.1); // planet
    struct reb_orbit o =  reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
    printf("a=%f e=%f\n", o.a, o.e);
    ```
    The last argument of the `reb_tools_particle_to_orbit` function is the primary particle, i.e. the star or the centre of mass. 
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1)  # star
    sim.add(a=1, e=0.1) # planet
    o = sim.particles[1].calculate_orbit(primary=sim.particles[0])
    print(o.a, o.e)
    ```
    If `primary` is not given, Jacobi coordinates are used.
    
    You can also calculate the orbits of all particles in the simulation. 
    ```python
    sim = rebound.Simulation()
    sim.add(m=1)  # star
    sim.add(a=1, e=0.1) # planet
    sim.add(a=2, e=0.1) # planet
    orbits = sim.calculate_orbits()
    for o in orbits:
        print(o.a, o.e)
    ```


## Conversion functions
### True anomaly

The following function returns the true anomaly $f$ for a given eccentricity $e$ and mean anomaly $M$:
=== "C"
    ```c
    double f = reb_tools_M_to_f(0.1, 1.); // e=0.1, M=1.0
    ```

=== "Python"
    ```python
    f = rebound.M_to_f(0.1, 1.0) # e=0.1, M=1.0
    ```


The following function returns the true anomaly $f$ for a given eccentricity $e$ and eccentric anomaly $E$:
=== "C"
    ```c
    double f = reb_tools_E_to_f(0.1, 1.); // e=0.1, E=1.0
    ```

=== "Python"
    ```python
    f = rebound.E_to_f(0.1, 1.0) # e=0.1, E=1.0
    ```

### Eccentric anomaly

The following function returns the eccentric anomaly $E$ for a given eccentricity $e$ and mean anomaly $M$:
=== "C"
    ```c
    double f = reb_tools_M_to_E(0.1, 1.); // e=0.1, M=1.0
    ```

=== "Python"
    ```python
    f = rebound.M_to_E(0.1, 1.0) # e=0.1, M=1.0
    ```
