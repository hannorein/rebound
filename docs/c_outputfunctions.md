# C output functions

The functions listed here provide various output functionality. 

## Output check

```c
int reb_simulation_output_check(struct reb_simulation* r, double interval);
```

This function can be used to trigger outputs at regular time intervals.
The function returns 1 if an output is required and 0 otherwise.
Typically, you would use this within the heartbeat function to generate equally spaced outputs as in this example:

```c
void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 100.)){
        printf("t = %f\n", r->t); // Will print current time every 100 time units.
    }
}
```

## Timing

```c
void reb_simulation_output_timing(struct reb_simulation* r, const double tmax);
```

This function outputs various status information on the screen. An example output looks as follows

```
N_tot=  10   t= 4610.00   dt= 4.000  cpu= 0.002472 [s]  t/tmax= 0.01%
```

It shows the number of particles, the current time and time-step, as well as the time since the last output. If `tmax` is non-zero, then the last number indicates how far the simulation has progressed. 

## ASCII orbits 
```c
void reb_simulation_output_orbits(struct reb_simulation* r, char* filename);
```
This function creates or appends an ASCII file with orbital parameters of all particles.
The orbital parameters are calculated in Jacobi coordinates.
Particles are assumed to be sorted from the inside out, the central object having index 0. 
Each time the function is called N-1 rows are appended to the file with name filename.
Each row in the file corresponds to one particle and contains the following columns (tab separated):

 - time
 - semi-major axis
 - eccentricity
 - inclination
 - Omega (longitude ascending node)
 - omega (argument of pericenter)
 - lambda (mean longitude)
 - orbital period,
 - f (true anomaly) 

## ASCII coordinates
```c
void reb_simulation_output_ascii(struct reb_simulation* r, char* filename);
```
This function creates or appends an ASCII file with the positions and velocities of all particles to an ASCII file.

## Velocity dispersion
```c
void reb_simulation_output_velocity_dispersion(struct reb_simulation* r, char* filename);
```
This function creates or appends an ASCII file with the current velocity dispersion of all particles. 
This is useful for ring simulations where one wants to monitor that the system has reached an equilibrium.

## Binary snapshot
```c
void reb_simulation_save_to_file(struct reb_simulation* r, const char* filename);
```

These functions save the `reb_simulation` structure as a binary file.
It can be used to save the current status of a REBOUND simulation and later restart the simulation.
If the file exists, this function will append a snapshot.

