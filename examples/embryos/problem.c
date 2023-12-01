/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "integrator_trace.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
int Nparticles = 30;
double tmax = 1e4 * 2 * M_PI;

char title[100] = "bs_embryo_col";
char title_remove[100] = "rm -rf bs_embryo_col";

int main(int argc, char* argv[]){
    system(title_remove);
    struct reb_simulation* r = reb_create_simulation();
//    r->integrator = REB_INTEGRATOR_BS;

    r->dt = (5./365.) * 2. *M_PI;
    r->softening = 3e-8;
    r->integrator = REB_INTEGRATOR_BS;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = reb_collision_resolve_merge;        // Choose merger collision routine.
    r->ri_ias15.adaptive_mode = 2;
    //r->ri_mercurius.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    //r->ri_trace.S_peri = reb_integrator_trace_switch_fdot_peri;
    //r->ri_trace.vfac_p = 16.0;

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
      reb_add_fmt(r, "primary m r a e Omega omega f", earth, m, 1.16e-5, a, e, Omega, omega, f);
    }

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    e_init = reb_tools_energy(r);
    printf("Initial Energy: %e\n", e_init);
    clock_t begin = clock();
    reb_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\nFinal Energy: %e\n", reb_tools_energy(r));
    printf("Conservation: %f\n", (reb_tools_energy(r) - e_init) / e_init);
    printf("Time Spent: %f\n", time_spent);
    reb_free_simulation(r);
}

void heartbeat(struct reb_simulation* r){
  if (reb_output_check(r, 10.*2.*M_PI)){
      reb_output_timing(r, 1e4 * 2 * M_PI);
  }

  if (reb_output_check(r, 1. * 2.*M_PI)){

    FILE* f = fopen(title, "a");
    fprintf(f, "%e,%e,%d\n", r->t, fabs((reb_tools_energy(r) - e_init) / e_init),r->N);
    fclose(f);
  }
}
