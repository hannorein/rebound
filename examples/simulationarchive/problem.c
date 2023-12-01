/**
 * Simulationarchive
 *
 * This example shows how to use the Simulationarchive.
 * We integrate a two planet system forward in time using
 * the WHFast integrator. The simulation can be interrupted
 * at any time. On the next run, the program will try to reload
 * the latest data from the Simulationarchive. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"


int main(int argc, char* argv[]) {
    char* filename = "simulationarchive.bin";

    // Trying to open a Simulationarchive file
    struct reb_simulationarchive* sa = reb_simulationarchive_create_from_file(filename);
    if (sa==NULL){
        printf("Can not open file.\n");
    }
    // Get a simulation from the file (if possible, otherwise NULL is returned)
    struct reb_simulation* r = reb_simulation_create_from_simulationarchive(sa,-1);
    // Whenever you've opened a Simulationarchive and don't need it anymore, close it.
    reb_simulationarchive_free(sa);
    // Check if we were successful
    if (r==NULL){
        printf("No simulationarchive found. Creating new simulation.\n");
        r= reb_simulation_create();
        reb_simulation_add_fmt(r, "m", 1.0);                   // star
        reb_simulation_add_fmt(r, "m a e", 1e-3, 1.0, 0.01);   // planet 1
        reb_simulation_add_fmt(r, "m a e", 1e-3, 2.3, 0.01);   // planet 2
        reb_simulation_move_to_com(r);
        r->dt = 6./365.25*2.*M_PI;              // 6 days in units where G=1 
        r->ri_whfast.safe_mode = 0;             // The Simulationarchive works with both safe_mode on and off           
        r->ri_whfast.corrector = 5;    
        r->integrator = REB_INTEGRATOR_WHFAST;
    }else{
        printf("Found simulationarchive. Loaded snapshot at t=%.16f.\n",r->t);
    }
    
    // Automatically create a snapshot every 100 time units
    reb_simulation_save_to_file_interval(r,filename,100.);
    // Alternatively, you can also create a snapshot every 5 seconds (walltime)
    //reb_simulation_save_to_file_walltime(r,filename,5.);

    // Run the integration (this will be very quick in this example)
    reb_simulation_integrate(r, r->t+2000); // integrate (a little further than where we currently are)
    printf("Final time: %f\n",r->t);
    
    // You can also manually append a snapshot
    reb_simulation_save_to_file(r,filename);

    // Free the simulation to free up memory
    reb_simulation_free(r);
}


