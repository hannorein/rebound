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
int Nparticles = 30;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
//    r->integrator = REB_INTEGRATOR_BS;

    r->dt = (5./365.) * 2. *M_PI;
    r->softening = 3e-8;
    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_tr.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    r->ri_tr.S_peri = reb_integrator_trace_switch_fdot_peri;
    r->ri_tr.vfac_p = 16.0;

    // Initialize masses
    struct reb_particle earth = {0};
    earth.m = 1;
    reb_add(r, earth);

    r->rand_seed = 1;
    r->heartbeat = heartbeat;
    // Test particles
    for (unsigned int i = 0; i < Nparticles; i++){
      double m = reb_random_uniform(r, 0.6*3.694e-8, 0.2*3e-6);
      double a = reb_random_uniform(r, 0.5, 1.2);
      double e = reb_random_uniform(r, 0.0, 0.01);
      double Omega = reb_random_uniform(r, 0.0, 2 * M_PI);
      double omega = reb_random_uniform(r, 0.0, 2 * M_PI);
      double f = reb_random_uniform(r, 0.0, 2 * M_PI);
      reb_add_fmt(r, "primary m a e Omega omega f", earth, m, a, e, Omega, omega, f);
    }

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    e_init = reb_tools_energy(r);
    printf("Initial Energy: %e\n", e_init);
    reb_integrate(r, 1e4 * 2 * M_PI);
    printf("\nFinal Energy: %e\n", reb_tools_energy(r));
    printf("Conservation: %f\n", (reb_tools_energy(r) - e_init) / e_init);
    reb_free_simulation(r);
}

void heartbeat(struct reb_simulation* r){
  if (reb_output_check(r, 10.*2.*M_PI)){
      reb_output_timing(r, 1e4 * 2 * M_PI);
  }
}
