#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <stdbool.h>
#include <gmp.h>

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

    reb_simulation_set_integrator(r, "whfast");
    //reb_simulation_set_integrator(r, "whfast512");
    //struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
    //whfast512->gr_potential = 0;
    return r;
}
extern void reb_integrator_whfast512_kepler_step(struct reb_simulation* const r, int N_steps);

   void init_mpfr();
int main(int argc, char* argv[]) {
    mpf_set_default_prec(200);
    init_mpfr();
    double a = 0.4;
    double e = 0.2;
    setbuf(stdout, NULL);
    struct reb_simulation* r = setup_sim(a,e);
    for (size_t i = 100; r->t<1e9*M_PI*2.0; i*=1.01){
        r->dt = 5.0/365.25*2*M_PI*(1.0+((double)i)/1e15);
        //reb_integrator_whfast512_kepler_step(r, i);
        for (size_t s=0; s<i; s++){
            reb_integrator_whfast_kepler_solver(&r->particles[1], 1.0, r->dt, r);
            r->t += r->dt;
        }
        struct reb_orbit o = reb_orbit_from_particle(1., r->particles[1], r->particles[0]);
        //printf("%e %e %e \n", r->t, fabs((a-o.a)/a), fabs((e-o.e)/e));
        //fflush(stdout);
        
    }
    return 1;
}
