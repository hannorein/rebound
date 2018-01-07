/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
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
    struct reb_simulation* r = reb_create_simulation();
    r->dt = 0.1;
    r->heartbeat = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.

    struct reb_particle p1 = {0}; // always initizialize a struct with this syntax to ensure all variables are set to 0.
    p1.m = 1.;
    reb_add(r, p1);  // reb_add makes a copy of the particle and adds it to the simulation.
    
    struct reb_particle p2 = {0};
    p2.x = 1;
    p2.vy = 1;
    reb_add(r, p2); // notice that we didn't set a mass. All parameters default to 0.

    reb_integrate(r,100.);
}

