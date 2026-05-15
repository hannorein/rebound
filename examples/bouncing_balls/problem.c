/**
 * Bouncing balls
 * 
 * This example is a simple test of collision detection 
 * methods.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
    
    // Setup constants
    reb_simulation_set_integrator(r, "leapfrog");
    r->gravity    = REB_GRAVITY_BASIC;
    r->collision    = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->usleep    = 1000;            // Slow down integration (for visualization only)
    r->dt = 1e-2;

    r->root_size = 3;
    
    // Initial conditions
    {
        struct reb_particle p = {0};
        p.x  = 1; p.y  = 1; p.z  = 1;
        p.m  = 1;
        p.r  = 0.1;
        reb_simulation_add(r, p);
    }
    {
        struct reb_particle p = {0};
        p.x  = -1; p.y  = -1; p.z  = -1;
        p.m  = 1;
        p.r  = 0.1;
        reb_simulation_add(r, p);
    }

    reb_simulation_integrate(r, INFINITY);
}

