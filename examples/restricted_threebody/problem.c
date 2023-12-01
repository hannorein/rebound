/**
 * Restricted three body problem.
 *
 * This example simulates a disk of test particles around 
 * a central object, being perturbed by a planet. 
 * It uses the heliocentric version of WHFast. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    // You can turn orbits on and off by pressing `w`.
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator   = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    r->boundary     = REB_BOUNDARY_OPEN;
    r->softening    = 1e-6;
    r->dt           = 1.0e-2*2.*M_PI;
    r->N_active     = 2;    // Only the star and the planet have non-zero mass
    r->heartbeat    = heartbeat;

    reb_simulation_configure_box(r,8.,1,1,1); // Box with size 8 AU
    
    // Initial conditions for star
    struct reb_particle star = {0};
    star.m  = 1;
    reb_simulation_add(r, star);

    // Initial conditions for planet
    double planet_e = 0.;
    struct reb_particle planet = {0};
    planet.x  = 1.-planet_e;
    planet.vy = sqrt(2./(1.-planet_e)-1.);
    planet.m  = 1e-2;
    reb_simulation_add(r, planet);
    reb_simulation_move_to_com(r);
    
    while(r->N<10000){
        double x    = ((double)rand()/(double)RAND_MAX-0.5)*8.;
        double y    = ((double)rand()/(double)RAND_MAX-0.5)*8.;
        double a    = sqrt(x*x+y*y);
        double phi  = atan2(y,x);
        if (a<.1) continue;
        if (a>4.) continue;

        double vkep = sqrt(r->G*star.m/a);
        struct reb_particle testparticle = {0};
        testparticle.x  = x;
        testparticle.y  = y; 
        testparticle.z  = 1.0e-2*x*((double)rand()/(double)RAND_MAX-0.5);
        testparticle.vx = -vkep*sin(phi);
        testparticle.vy = vkep*cos(phi);
        reb_simulation_add(r, testparticle);
    }

    reb_simulation_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 20.*M_PI)){
        reb_simulation_output_timing(r, 0);
        reb_simulation_output_orbits(r, "orbit.txt");
    }
}
