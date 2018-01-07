/**
 * A string of solid spheres bouncing
 * 
 * This example tests collision detection methods.
 * The example uses a non-square, rectangular box. 10 particles are placed
 * along a line. All except one of the particles are at rest initially.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

int main(int argc, char* argv[]){
    struct reb_simulation* const r = reb_create_simulation();
    // Setup constants
    r->dt             = 1e-3;
    r->integrator        = REB_INTEGRATOR_LEAPFROG;
    r->boundary        = REB_BOUNDARY_PERIODIC;
    r->collision        = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->gravity        = REB_GRAVITY_NONE;
    r->usleep        = 5000;            // Slow down integration (for visualization only)
    
    reb_configure_box(r,10.,3,1,1);  // boxsize 10., three root boxes in x direction, one in y and z
    r->nghostx = 1; 
    r->nghosty = 1; 
    r->nghostz = 0;

    // Initial conditions
    for(int i=0;i<10;i++){
        struct reb_particle p = {0};
        p.x  = -r->boxsize.x/2.+r->boxsize.x*(double)i/10.; p.y  = 0; p.z  = 0;
        p.m  = 1;
        p.r  = 1;
        reb_add(r, p);
    }

    // Give one particle a kick
    r->particles[0].vx = 20;

    reb_integrate(r,INFINITY);
}
