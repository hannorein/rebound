/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "integrator_trace.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
int Nparticles = 999;
// accretion of the moon
int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
//    r->integrator = REB_INTEGRATOR_BS;

    r->dt = 0.2;
    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_tr.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    r->ri_tr.S_peri = reb_integrator_trace_switch_fdot_peri;
    r->ri_tr.vfac_p = 16.0;

    // Initialize masses
    struct reb_particle earth = {0};
    earth.m = 1;
    reb_add(r, earth);

    r->rand_seed = 1;
    // Test particles
    for (unsigned int i = 0; i < Nparticles; i++){
      double m = reb_random_uniform(r, 3.2e-7, 3.2e-4);
      double a = reb_random_uniform(r, 0.7, 1.3);
      double inc = reb_random_uniform(r, 0.0, 50. * M_PI/180.);
      reb_add_fmt(r, "primary m a inc", earth, m, a, inc);
    }

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    e_init = reb_tools_energy(r);
    printf("Initial Energy: %e\n", e_init);
    reb_integrate(r, 6 * M_PI);
    printf("Final Energy: %e\n", reb_tools_energy(r));
    printf("Conservation: %f\n", (reb_tools_energy(r) - e_init) / e_init);
    reb_free_simulation(r);
}

void heartbeat(struct reb_simulation* r){
}
