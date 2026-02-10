/**
 * REBOUND API without simulations
 * 
 * REBOUND is a shared library. This means you can not only run full simulations
 * from your program, but also use internal functions from REBOUND in your
 * own programs. Here we demonstrate some common use cases: We use a Kepler
 * solver to move a planet along it's orbit. And we calculate orbital elements
 * for a pair of particles. Note that we never initialize a simulation.
 *
 *
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    // Let's use the WHFAST's built-in Kepler solver to move a single particle
    // along its Keplerian orbit. This solver is very fast and accurate for 
    // eccentric as well as hyperbolic orbits.

    // We define a massless particle on a circular orbit:
    struct reb_particle p = {.x=1,.vy=1};
    printf("Initial position: %f %f %f\n", p.x, p.y, p.z);

    // Now we move the particle forward along it's orbit.
    // The central object is at the coordinate origin and has a mass of 1.0.
    double G = 1.0;     // Working in units where G=1, but you can choose other units if you prefer.
    double GM = G*1.0;  // The gravitational parameter (G*mass);
    double dt = M_PI;   // Time interval, here: half an orbital period
    reb_whfast_kepler_solver(NULL, &p, GM, 0,  dt); 
    // We pass NULL to indicate that particle p is not part of a REBOUND simulation
    // The second argument can be a pointer to a list of particles. Here we just use the 0th entry. 

    // The particle is now on the opposite site of the primary.
    printf("Final position: %f %f %f\n", p.x, p.y, p.z);    

    
    // Let's use the built-in coordinate transformations tools in REBOUND
    // to calculate a particle's orbital parameters. These routines are well
    // tested and work in all edge cases.

    // We define the primary with respect to which we will calculate the orbital
    // parameters. This could be an actual particle (e.g. the Sun in the Solar System)
    // or a virtual partical, for example when calculating Jacobi coordinates.
    struct reb_particle primary = {.m=1}; // Particle with mass 1, at origin, at rest
    int err = 0; // Error code.
    struct reb_orbit o = reb_orbit_from_particle_err(G, p, primary, &err);

    // Checking if an error occured.
    // This can happen if the primary has no mass or if the particles are on top of each other. 
    if (err){
        printf("An error occured during orbit calculation.\n");    
    }

    // Printing out some parameters (see rebound.h for all parameters)
    printf("Semi-major axis: %f\nEccentricity: %f\nInclination: %f\n", o.a, o.e, o.inc);    
   
    // We can also for example calculate the eccentric anomaly for a given eccentricity and mean anomaly
    printf("Eccentric anomaly: %f\n", reb_M_to_E(o.e, o.M));    
   
}

