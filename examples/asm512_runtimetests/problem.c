/**
 * Unit tests for WHFast512
 *
 * This file contains units tests for WHFast512.
 * Note that these are not run automatically 
 * because GitHub's CI does not support AVX5212.
 */

#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sched.h>
#include <stdbool.h>
#include <sys/time.h>

struct reb_simulation* setup_sim(char* integrator){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 6.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
   
    reb_simulation_add_fmt(r, "solarsystem");

    reb_simulation_set_integrator(r, integrator);
    if (strcmp(integrator,"whfast512")==0){
        struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
        whfast512->gr_potential = 0;
        whfast512->concatenate_steps = 1e6;
    }
    if (strcmp(integrator,"whfast")==0){
        struct reb_integrator_whfast_state* whfast = r->integrator.state;
        whfast->safe_mode = 0;
    }
    return r;
}

int main(int argc, char* argv[]) {
    struct timeval time_beginning;
    struct timeval time_end;
    double tmax_years = 1e5;
    double walltime[3];

    for (int i=0; i<2; i++){
        char* integrator = "whfast512";
        if (i==1) integrator = "whfast";
        printf("###################################\n");
        printf("integrator:     %s\n", integrator);
        struct reb_simulation* r = setup_sim(integrator);
        gettimeofday(&time_beginning,NULL);
        reb_simulation_integrate(r, tmax_years*M_PI*2.0);
        gettimeofday(&time_end,NULL);
        walltime[i] = (time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6)/(r->t/2/M_PI) * 5e9/60./60.;
        printf("time to 5 Gyr:  %.8f hours\n", walltime[i]);
        
        reb_simulation_free(r);
    }

    printf("###################################\n");
    printf("speedup asm512/whfast    = %.5fx\n", walltime[1]/walltime[0]);

    return 1;
}
