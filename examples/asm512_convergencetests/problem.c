#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <stdbool.h>
#include <string.h>

void setup_sim(char* integrator, int corrector, double dt){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = dt/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
    reb_simulation_add_fmt(r, "solarsystem");
    reb_simulation_set_integrator(r, integrator);
    double tmax = 1e2*M_PI*2;
    long steps = ceil(tmax/r->dt);
    if (strcmp(integrator,"asm512")==0){
        struct reb_integrator_asm512_state* asm512 = r->integrator.state;
        asm512->gr_potential = 0;
        asm512->concatenate_steps = steps;
        asm512->corrector = corrector;
    }
    if (strcmp(integrator,"whfast512")==0){
        struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
        whfast512->gr_potential = 0;
    }
    if (strcmp(integrator,"whfast")==0){
        struct reb_integrator_whfast_state* whfast = r->integrator.state;
        whfast->safe_mode = 0;
        whfast->corrector = corrector;
    }
    double E0 = reb_simulation_energy(r);
    reb_simulation_integrate(r, tmax);
    double E1 = reb_simulation_energy(r);
    printf("%.16e ", fabs((E0-E1)/E0));
    reb_simulation_free(r);
}



extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);

int main(int argc, char* argv[]) {
    for (double dt = 1e-1; dt < 1e2; dt*=1.02){
        printf("%.16e ", dt);
        setup_sim("asm512", 0, dt);
        setup_sim("asm512", 17, dt);
        setup_sim("whfast", 0, dt);
        setup_sim("whfast", 17, dt);
        setup_sim("whfast512", 0, dt);
        printf("\n");
    }

}
