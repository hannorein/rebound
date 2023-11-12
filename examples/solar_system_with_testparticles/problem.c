/**
 * Solar System with test particles
 *
 * This example integrates all planets of the Solar
 * System and 10000 test particles. The initial data comes 
 * from the NASA HORIZONS system and was saved to
 * a binary file beforehand. The integrator used is WHFast
 * with a 4 day timestep. Note that close encounters are
 * not resolved. The OpenMP speedup you get depends on the 
 * compiler and CPU that you are using. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create_from_file("ss-2023-11-12.bin", 0);
    
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->dt           = 4./365.25*2.*M_PI;        // 4days
    r->integrator   = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates   = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    r->heartbeat    = heartbeat;
    r->N_active     = r->N;

    // Test Particles
    for (int i=0;i<10000;i++){
        double a = reb_random_uniform(r, 0.4,20.);
        double e = reb_random_uniform(r, 0.01,0.2);
        double omega = reb_random_uniform(r, 0.,2.*M_PI);
        double f = reb_random_uniform(r, 0.,2.*M_PI);
        struct reb_particle p = reb_particle_from_orbit(1.,r->particles[0],0.,a,e,0.,0.,omega,f);
        reb_simulation_add(r, p); 
    }
    reb_simulation_move_to_com(r);

    // Integrate forever
    reb_simulation_integrate(r, INFINITY);

    // cleanup
    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 100.)){
        reb_simulation_output_timing(r, INFINITY);
    }
}

