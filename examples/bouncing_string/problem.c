/**
 * A string of solid spheres bouncing
 * 
 * This example tests collision detection methods.
 * The example uses a non-square, rectangular box. 10 particles are placed
 * along a line. All except one of the particles are at rest initially.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    // Setup modules and constants
    r->dt                = 1e-3;
    reb_simulation_set_integrator(r,"leapfrog");
    r->boundary          = REB_BOUNDARY_PERIODIC;
    r->collision         = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->gravity           = REB_GRAVITY_NONE;
    r->usleep            = 5000;      // Slow down integration (for visualization only)
    
    // boxsize 10., three root boxes in x direction, one in y and z
    r->root_size = 10.0;
    r->N_root_x = 3; 

    r->N_ghost_x = 1; 
    r->N_ghost_y = 1; 
    r->N_ghost_z = 0;

    // Initial conditions
    struct reb_vec3d boxsize = {
        .x = r->root_size*(double)r->N_root_x,
        .y = r->root_size*(double)r->N_root_y,
        .z = r->root_size*(double)r->N_root_z,
    };
    for(int i=0;i<10;i++){
        struct reb_particle p = {0};
        p.x  = -boxsize.x/2.+boxsize.x*(double)i/10.; p.y  = 0; p.z  = 0;
        p.m  = 1;
        p.r  = 1;
        reb_simulation_add(r, p);
    }

    // Give one particle a kick
    r->particles[0].vx = 20;

    reb_simulation_integrate(r,INFINITY);
}
