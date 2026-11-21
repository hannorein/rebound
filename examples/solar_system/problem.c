/**
 * Solar System
 *
 * This example integrates all planets of the Solar
 * System. The data comes from the NASA HORIZONS system. 
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void heartbeat(struct reb_simulation* r);
double e_init;
double tmax;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->dt                   = 0.05*M_PI*2.0;       // in units of years/(2*pi)
    tmax                    = 1e6*M_PI*2.0;        // in units of years/(2*pi)
    reb_simulation_set_integrator(r, "whfast");
    // reb_simulation_set_integrator(r, "whfast");    // Alternative non-symplectic integrator
    struct reb_integrator_whfast_state* whfast = r->integrator.state;
    whfast->safe_mode  = 0;             // Turn off safe mode. Need to call reb_simulation_synchronize() before outputs. 
    whfast->corrector  = 11;            // 11th order symplectic corrector
    r->heartbeat            = heartbeat;
    r->exact_finish_time    = 1;     // Finish exactly at tmax in reb_simulation_integrate(). Default is already 1.

    // Initial conditions
    reb_simulation_add_fmt(r, "solar system"); // Built in dataset for testing.
    reb_simulation_move_to_com(r);

    e_init = reb_simulation_energy(r);
    remove("energy.txt");
    reb_simulation_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10*M_PI*2.0)){ // in units of years/(2*pi)
        reb_simulation_output_timing(r, tmax);
        reb_simulation_synchronize(r);
        FILE* f = fopen("energy.txt","ab");
        double e = reb_simulation_energy(r);
        fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        fclose(f);
    }
}

