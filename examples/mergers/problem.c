/**
 * Colliding and merging planets
 * 
 * This example integrates a densely packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. The IAS15 
 * integrator with adaptive timestepping is used. The bodies have a finite
 * size and merge if they collide. Note that the size is unphysically large
 * in this example. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    r->dt                   = 0.01*2.*M_PI;                // initial timestep
    r->integrator           = REB_INTEGRATOR_IAS15;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = reb_collision_resolve_merge; // Choose merger collision routine.
    r->heartbeat            = heartbeat;

    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.1;
    reb_simulation_add(r, star);
    
    // Add planets
    int N_planets = 7;
    for (int i=0;i<N_planets;i++){
        double a = 1.+(double)i/(double)(N_planets-1);   // semi major axis in AU
        double v = sqrt(1./a);                           // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 1e-4; 
        planet.r = 4e-2;                                 // radius in AU (it is unphysically large in this example)
        planet.last_collision = 0;                       // The first time particles can collide with each other
        planet.x = a; 
        planet.vy = v;
        reb_simulation_add(r, planet); 
    }
    reb_simulation_move_to_com(r);                       // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    reb_simulation_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.*2.*M_PI)){  
        reb_simulation_output_timing(r, 0);
    }
}
