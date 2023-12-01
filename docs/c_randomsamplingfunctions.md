# Random sampling

REBOUND includes several functions which help sampling random numbers from various distributions.
Each function takes a pointer to a `struct reb_simulation` as a first argument.
This is because the simulation structure stores the random number generator state in the member `rand_seed`. 
When a simulation is created, `rand_seed` is initialized using the current time and process id. 
These functions are thread safe.
The `rand_seed` variable is stored in binary files. This makes the random number generator reproducible which can be very helpful when debugging simulations that use random numbers. 

The following example draws a number in the interval between 0 and $2\pi$ from a uniform distribution.
```c
struct reb_simulation* r = reb_simulation_create();
double phi = reb_random_uniform(r, 0., 2.*M_PI);
```

# Uniform
This function returns a uniformly distributed random variable between `min` and `max`.
```c
double reb_random_uniform(struct reb_simulation* r, double min, double max);
```

# Power law
This function returns a random variable drawn form a power law distribution with slope `slope` between `min` and `max`. 
```c
double reb_random_powerlaw(struct reb_simulation* r, double min, double max, double slope);
```

# Normal
This function returns a random number drawn from a normal distribution centerd on zero and with variance `variance`.
It uses the algorithm by D.E. Knut, 1997, The Art of Computer Programming, Addison-Wesley. 
```c
double reb_random_normal(struct reb_simulation* r, double variance);
```

# Rayleigh
This function returns a random variable drawn form a Rayleigh distribution with scale parameter `sigma`. 
```c
double reb_random_rayleigh(struct reb_simulation* r, double sigma);
```


