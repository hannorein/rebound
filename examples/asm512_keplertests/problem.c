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
    reb_simulation_add_fmt(r, "a e ", a, e);
    //reb_simulation_add_fmt(r, "a e uniform(f)", a, e);
    for (int i=0; i<7; i++){
        reb_simulation_add_fmt(r, "a", a+1+i);
    }

    reb_simulation_set_integrator(r, "asm512");
    struct reb_integrator_asm512_state* asm512 = r->integrator.state;
    asm512->gr_potential = 0;
    return r;
}

extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);


int main(int argc, char* argv[]) {
    int Na = 100;
    int Ne = 100;
    for (int ia=0; ia<Na; ia++){
        for (int ie=0; ie<Ne; ie++){
            double a = 0.1 + (0.95-0.1)*(double)ia/(double)(Na+1);
            double e = 0.0 + (0.95-0.0)*(double)ie/(double)(Ne+1);
            struct reb_simulation* r = setup_sim(a,e);
            reb_simulation_set_integrator(r, "asm512");
            reb_integrator_asm512_kepler_step(r, 730); // 10yrs
            struct reb_orbit o = reb_orbit_from_particle(1., r->particles[1], r->particles[0]);
            struct reb_particle p;
            //p = r->particles[0];
            //printf("[0] x y z = %f %f %f   %f %f %f\n", p.x, p.y, p.z, p.vx, p.vy, p.vz);
            //p = r->particles[1];
            //printf("[1] x y z = %f %f %f   %f %f %f\n", p.x, p.y, p.z, p.vx, p.vy, p.vz);
            printf("%e %e %e %e\n", a, e, fabs((a-o.a)/a), a-o.a);
            reb_simulation_free(r);
        }
    }
    return 1;
}
