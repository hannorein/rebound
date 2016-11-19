/**
 * Simulation Archive
 *
 * This example shows how to use the Simulation Archive.
 * We integrate a two planet system forward in time using
 * the WHFast integrator. The simulation can be interrupted
 * at any time. On the next run, the program will try to reload
 * the latest data from the Simulation Archive. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "rebound.h"


int main(int argc, char* argv[]) {
    char filename[512] = "simulationarchive.bin";
    // Trying to restart from the Simulation Archive.
    struct reb_simulation* r = reb_simulationarchive_restart(filename);
    // Check if that was successful
    if (r==NULL){
        printf("No simulation archive found. Creating new simulation.\n");
        r= reb_create_simulation();
        struct reb_particle star = {.m=1.};
        reb_add(r, star);
        struct reb_particle planet1 = reb_tools_orbit2d_to_particle(r->G, star, 1e-3, 1., 0.01, 0., 0.);
        reb_add(r, planet1);
        struct reb_particle planet2 = reb_tools_orbit2d_to_particle(r->G, star, 1e-3, 2.3, 0.01, 0., 0.);
        reb_add(r, planet2);
        reb_move_to_com(r);
        r->dt = 6./365.25*2.*M_PI;                      // 6 days in units where G=1 
        r->simulationarchive_interval = 2.*M_PI*100000.; // output data every 100000 years
        r->ri_whfast.safe_mode = 0;                      
        r->ri_whfast.corrector = 5;    
        r->integrator = REB_INTEGRATOR_WHFAST;
    }else{
        printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
    }

    r->simulationarchive_filename = filename;

    reb_integrate(r, 2.*M_PI*1e10); //10 Gyr
}


