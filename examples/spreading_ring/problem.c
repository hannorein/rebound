/**
 * Spreading ring
 *
 * A narrow ring of collisional particles is spreading.
 * An error message will alert you to the fact that we
 * do not take the mass of particles into account when
 * calculating mutual gravity. However, we use the mass
 * during the collision resolve phase. Thus you can ignore 
 * the warning in this case.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator        = REB_INTEGRATOR_LEAPFROG;
    r->collision         = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->boundary          = REB_BOUNDARY_OPEN;
    r->G                 = 1;        
    r->N_active          = 1;
    r->softening         = 0.01;        
    r->dt                = 1e-3;
    r->heartbeat         = heartbeat;
    
    double boxsize = 4.8;
    reb_simulation_configure_box(r, boxsize, 1, 1, 1);

    // Setup particles
    int _N = 1000;
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
    star.r        = 0.01;
    reb_simulation_add(r, star);

    while(r->N<_N){
        struct reb_particle pt = {0};
        double a     = reb_random_powerlaw(r, boxsize/2.9,boxsize/3.1,.5);
        double phi   = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.0001);
        double vkep  = sqrt(r->G*star.m/a);
        pt.vx        =  vkep * sin(phi);
        pt.vy        = -vkep * cos(phi);
        pt.m         = 0.0001;
        pt.r         = .3/sqrt((double)_N);
        reb_simulation_add(r, pt);
    }

    reb_simulation_integrate(r, INFINITY);

    // Cleanup
    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 0.0*r->dt)){
        reb_simulation_output_timing(r, 0);
    }
}
