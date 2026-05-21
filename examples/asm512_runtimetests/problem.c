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
    if (strcmp(integrator,"asm512")==0){
        struct reb_integrator_asm512_state* asm512 = r->integrator.state;
        asm512->gr_potential = 0;
        asm512->concatenate_steps = 1e6;
    }
    if (strcmp(integrator,"whfast512")==0){
        struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
        whfast512->gr_potential = 0;
    }
    if (strcmp(integrator,"whfast")==0){
        struct reb_integrator_whfast_state* whfast = r->integrator.state;
        whfast->safe_mode = 0;
    }
    return r;
}

extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);


int main(int argc, char* argv[]) {
    struct timeval time_beginning;
    struct timeval time_end;
    double tmax_years = 1e5;

    for (int i=0; i<3; i++){
        char* integrator = "asm512";
        if (i==1) integrator = "whfast512";
        if (i==2) integrator = "whfast";
        printf("###################################\n");
        printf("integrator:     %s\n", integrator);
        struct reb_simulation* r = setup_sim(integrator);
        gettimeofday(&time_beginning,NULL);
        reb_simulation_integrate(r, tmax_years*M_PI*2.0);
        gettimeofday(&time_end,NULL);
        double walltime = time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
        printf("time:           %.8f seconds\n", walltime);
        printf("time to 5 Gyr:  %.8f hours\n", walltime/(r->t/2/M_PI) * 5e9/60./60.);
        
        reb_simulation_free(r);
    }
    return 1;
}
