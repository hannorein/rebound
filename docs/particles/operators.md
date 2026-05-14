# Operators

## Adding, subtracting, multiplying particles
REBOUND allows you to multiply particles with scalars.
In the code blow, the particle's position and velocity coordinates, and its mass are all multiplied by a constant.

=== "C"
    ```c
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_add_fmt(r, "m x vx", 1., 1., 1.);
    reb_particle_imul(&(r->particles[0]), 2.);
    ```
=== "Python"
    ```python
    sim = rebound.Simulation()
    sim.add(m=1., x=1., vx=1.)
    sim.particles[0] *= 2.
    ```
You can also add or subtract particles from each other. 
This will add or subtract the particles' position, velocity and mass from each other.

=== "C"
    ```c
    struct reb_particle p1 = {.m=1., .x=1, .vx=1};
    struct reb_particle p2 = {.m=2., .x=2, .vx=2};
    reb_particle_iadd(&p1, &p2); // p1.m, p1.x, p1.vx will now all be 3. p2 remains unchanged.
    reb_particle_isub(&p1, &p2); // p1.m, p1.x, p1.vx will be 1.
    ```
=== "Python"
    ```python
    p1 = rebound.Particle(m=0., x=1., vx=1.)
    p2 = rebound.Particle(m=0., x=1., vx=1.)
    p1 += p2 # p1.m, p1.x, p1.vx will now all be 3. p2 remains unchanged.
    p1 -= p2 # p1.m, p1.x, p1.vx will be 1.
    ```

In all the C functions, the first particle gets modified in place. 
In python, one can also use the multiply, add, and subtract operations to create new particles.
The following operations do not affect the original particles `p1` and `p2`.

```python
p1 = rebound.Particle(m=0., x=1., vx=1.)
p2 = rebound.Particle(m=0., x=1., vx=1.)
p3 = p1 + p2  # p3 is a new particle
p4 = p1 + p2  # p4 is a new particle
p5 = 2.*p1    # p5 is a new particle
```

!!! Tip
    These particle operations can be very helpful when initializing particles.
    For example, you can create initialize two particles using orbital parameters and then easily create another particle exactly in between these two particles.
    ```python
    sim = rebound.Simulation()
    sim.add(m=1)        # star
    sim.add(a=1)        # planet 1
    sim.add(a=2, f=0.1) # planet 2
    p_middle = (p1+p2)/2. # exactly in between the two planets
    ```

The distance between two particles in 3D space is often required in various calculations.
REBOUND has a convenience function for that:

=== "C"
    ```c
    struct reb_particle p1 = {.m=1., .x=1, .vx=1};
    struct reb_particle p2 = {.m=2., .x=2, .vx=2};
    double distance = reb_particle_distance(&p1, &p2); 
    ```
=== "Python"
    In python this is implemented using the power operator (`**`):
    ```python
    p1 = rebound.Particle(m=0., x=1., vx=1.)
    p2 = rebound.Particle(m=0., x=1., vx=1.)
    distance = p1 ** p2
    ```

