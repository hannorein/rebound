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

struct reb_simulation* setup_sim(double a, double e){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 5.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
   
    reb_simulation_add_fmt(r, "m", 1.0);
    for (int i=0; i<8; i++){
        reb_simulation_add_fmt(r, "a e uniform(f)", a, e);
    }

    reb_simulation_set_integrator(r, "asm512");
    struct reb_integrator_asm512_state* asm512 = r->integrator.state;
    asm512->gr_potential = 0;
    return r;
}

extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);


int main(int argc, char* argv[]) {
    int Na = 100;
    int Ne = 100;
    for (int ia=0; ia<Na; ia++){
        for (int ie=0; ie<Ne; ie++){
            double a = 0.05 + (0.95-0.05)*(double)ia/(double)(Na-1);
            double e = 0.0 + (0.99-0.0)*(double)ie/(double)(Ne-1);
            struct reb_simulation* r = setup_sim(a,e);
            reb_simulation_set_integrator(r, "asm512");
            int Nsteps = 1; // 10 years
            reb_integrator_asm512_kepler_step(r, Nsteps);
            int test_p = rand_r(&r->rand_seed) % 8;
            struct reb_orbit o = reb_orbit_from_particle(1., r->particles[test_p+1], r->particles[0]);
            double s = 0;
            if (a-o.a>0.0){
                s = 1.0;
            }
            if (a-o.a<0.0){
                s = -1.0;
            }
            uint64_t counter  = reb_asm512_counter(r, test_p);
            printf("%e %e %e %e %e\n", a, e, fabs((a-o.a)/a), s, ((double)counter)/((double)Nsteps));
            reb_simulation_free(r);
        }
    }
    return 1;
}
