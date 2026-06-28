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

struct reb_simulation* setup_sim(double a, double e, double de){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 5.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
   
    reb_simulation_add_fmt(r, "m", 1.0);
    for (int i=0; i<8; i++){
        reb_simulation_add_fmt(r, "a e uniform(f)", a, e+de*i);
    }

    reb_simulation_set_integrator(r, "whfast512");
    struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    whfast512->gr_potential = 0;
    return r;
}

extern void reb_integrator_whfast512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_whfast512_counter(struct reb_simulation* r, int test_p);


int main(int argc, char* argv[]) {
    int Na = 8*12;
    int Ne = 8*12;
    int id = -1;
    int ia_start = 0;
    int ia_end = Na;
    if (argc>=2){
        id = atoi(argv[1]);
        ia_start = id;
        ia_end = id+1;
    }
    char filename[1024];
    sprintf(filename, "out_parallel_%2d.txt", id);
    FILE* f = fopen(filename, "w");
    for (int ia=ia_start; ia<ia_end; ia++){
        for (int ie=0; ie<Ne; ie+=8){
            double a = 0.05 + (0.95-0.05)*(double)ia/(double)(Na-1);
            double de = (0.99-0.0)/(double)(Ne-1);
            double e = 0.0 + de*(double)ie;
            int Nsteps = 1e4*365.25/5.0;
            struct reb_simulation* r = setup_sim(a,e, de);
            reb_integrator_whfast512_kepler_step(r, Nsteps);
            for (int k=0; k<8; k++){
                struct reb_orbit o = reb_orbit_from_particle(1., r->particles[k+1], r->particles[0]);
                double s = 0;
                if (a-o.a>0.0){
                    s = 1.0;
                }
                if (a-o.a<0.0){
                    s = -1.0;
                }
                fprintf(f, "%e %e %e %e \n", a, e+k*de, fabs((a-o.a)/a), s);

            }
            reb_simulation_free(r);
        }
    }
    fclose(f);
    return 1;
}
