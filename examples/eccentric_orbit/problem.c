/**
 * Highly eccentric orbits
 *
 * This example uses the IAS15 integrator to simulate
 * a very eccentric planetary orbit. The integrator
 * automatically adjusts the timestep so that the pericenter passages
 * are resolved with high accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double timescale; // orbital timescale
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    // Setup constants
    r->G            = 1;        // Gravitational constant
    r->integrator   = REB_INTEGRATOR_IAS15;
    r->heartbeat    = heartbeat;
    r->ri_ias15.adaptive_mode = 2;    // Improved timestep criterion

    double e_testparticle = 1.-1e-7;    
    double mass_scale     = 1.;       // Some integrators have problems when changing the mass scale, IAS15 does not. 
    double size_scale     = 1.;       // Some integrators have problems when changing the size scale, IAS15 does not.

    struct reb_particle star = {0}; 
    star.m  = mass_scale;
    reb_simulation_add(r, star); 
    
    struct reb_particle planet; 
    planet.m  = 0;
    planet.x  = size_scale*(1.-e_testparticle); 
    planet.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle)*mass_scale/size_scale);
    reb_simulation_add(r, planet); 
    
    reb_simulation_move_to_com(r);
    
    // initial timestep
    r->dt     = 1e-13*sqrt(size_scale*size_scale*size_scale/mass_scale); 
    // calculate orbital timescale
    timescale = 2.*M_PI*sqrt(size_scale*size_scale*size_scale/mass_scale);

    reb_simulation_integrate(r, INFINITY);    
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r,timescale/0.1)){        // outputs to the screen every 0.1 orbits
        reb_simulation_output_timing(r, 0);
    }
}

