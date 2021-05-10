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
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation_from_binary("ss-2020-03-03.bin");
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
        struct reb_particle p = reb_tools_orbit_to_particle(1.,r->particles[0],0.,a,e,0.,0.,omega,f);
        reb_add(r, p); 
    }
    reb_move_to_com(r);
    reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 100.)){
        reb_output_timing(r, INFINITY);
    }
}

