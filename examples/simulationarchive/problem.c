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
    char* filename = "simulationarchive.bin";

    // Trying to open a SimulationArchive file
    struct reb_simulationarchive* sa = reb_open_simulationarchive(filename);
    if (sa==NULL){
        printf("Can not open file.\n");
    }
    // Get a simulation from the file (if possible, otherwise NULL is returned)
    struct reb_simulation* r = reb_create_simulation_from_simulationarchive(sa,-1);
    // Whenever you've opened a SimulationArchive and don't need it anymore, close it.
    reb_close_simulationarchive(sa);
    // Check if we were successful
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
        r->dt = 6./365.25*2.*M_PI;              // 6 days in units where G=1 
        r->ri_whfast.safe_mode = 0;             // The SimulationArchive works with both safe_mode on and off           
        r->ri_whfast.corrector = 5;    
        r->integrator = REB_INTEGRATOR_WHFAST;
    }else{
        printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
    }
    
    // Automatically create a snapshot every 100 time units
    reb_simulationarchive_automate_interval(r,filename,100.);
    // Alternatively, you can also create a snapshot every 5 seconds (walltime)
    //reb_simulationarchive_automate_walltime(r,filename,5.);

    // Run the integration (this will be very quick in this example)
    reb_integrate(r, r->t+2000); // integrate (a little further than where we currently are)
    printf("Final time: %f\n",r->t);
    
    // You can also manually append a snapshot
    reb_simulationarchive_snapshot(r,filename);

    // Free the simulation to free up memory
    reb_free_simulation(r);
}


