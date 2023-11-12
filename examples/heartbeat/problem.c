/**
 * How to use a heartbeat function
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units. We output the current time at every timestep.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

void heartbeat(struct reb_simulation* r){
    // This function gets called after every timestep.
    // Here, we simply print out the current simulation time.
    printf("%f\n",r->t);
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->dt = 0.1;
    r->heartbeat = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_simulation_integrate(). Default is already 1.

    struct reb_particle p1 = {0}; // always initizialize a struct with this syntax to ensure all variables are set to 0.
    p1.m = 1.;
    reb_simulation_add(r, p1);  // reb_simulation_add makes a copy of the particle and adds it to the simulation.
    
    reb_simulation_add_fmt(r, "a e", 1., 0.); // We can also add a particle using the reb_simulation_add_fmt function.

    reb_simulation_integrate(r,100.);

    reb_simulation_free(r);
}

