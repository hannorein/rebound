/**
 * Close Encounter
 *
 * This example integrates a densely packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. The IAS15 
 * integrator with adaptive timestepping is used. This integrator 
 * automatically decreases the timestep whenever a close 
 * encounter happens. IAS15 is very high order and ideally suited for the 
 * detection of these kind of encounters.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    r->dt         = 0.01*2.*M_PI;        // initial timestep
    r->integrator     = REB_INTEGRATOR_IAS15;
    r->heartbeat      = heartbeat;
    r->usleep    = 10000;        // Slow down integration (for visualization only)

    // Add star
    struct reb_particle star = {0};
    star.m = 1;
    reb_add(r, star);
    
    // Add planets
    int N_planets = 7;
    for (int i=0;i<N_planets;i++){
        double a = 1.+(double)i/(double)(N_planets-1);        // semi major axis
        double v = sqrt(1./a);                     // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 1e-4; 
        planet.x = a; 
        planet.vy = v; 
        reb_add(r, planet); 
    }
    reb_move_to_com(r);        // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10.*2.*M_PI)){  
        reb_output_timing(r, 0);
    }
}
