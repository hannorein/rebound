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
#include <sched.h>
#include <stdbool.h>

struct reb_simulation* setup_sim(){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 5.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
   
    reb_simulation_add_fmt(r, "solarsystem");

    return r;
}

extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);


int main(int argc, char* argv[]) {
    struct reb_simulation* r_asm = setup_sim();
    double E0 = reb_simulation_energy(r_asm);
    reb_simulation_set_integrator(r_asm, "asm512");
    struct reb_integrator_asm512_state* asm512 = r_asm->integrator.state;
    asm512->gr_potential = 0;
    asm512->concatenate_steps = 1e6;
    asm512->corrector = 17;
    for (double dT=1e1; r_asm->t<5e9*2*M_PI; dT=dT*1.05){
        reb_simulation_integrate(r_asm, r_asm->t+dT);
        double E1_asm = reb_simulation_energy(r_asm);
        char* mode = "a";
        if (dT==1e1) mode = "w";
        FILE* f = fopen("out.txt",mode);
        fprintf(f,"%e %e\n", r_asm->t, fabs((E0-E1_asm)/E0));
        fclose(f);
    }
    reb_simulation_free(r_asm);
    return 1;
}
