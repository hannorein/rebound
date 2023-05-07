/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;

int main(int argc, char* argv[]){

    // Initialize masses

    // Jupiter
    double jm = 9.55e-4;
    double ja = 5.2;
    double je = 0.05;

    // Saturn
    double sm = 2.857e-4;
    double sa = 9.58;
    double se = 0.95;
    double si = M_PI / 2.;

    int num_peris = 10;
    double peris[10] = {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 0.9, 2.};

    system("rm -rf energy_trace_grid.txt");
    FILE* f = fopen("energy_trace_grid.txt","a");
    fprintf(f, "p,e\n");
    for (int i = 0; i < num_peris; i++){
      err = 0.0;

      struct reb_simulation* r = reb_create_simulation();
      struct reb_particle star = {0};
      star.m = 1.;
      reb_add(r, star);
      reb_add_fmt(r, "m a e", jm, ja, je);

      se = (sa - peris[i]) / sa;
      printf("%e\n", se);
      reb_add_fmt(r, "m a e inc", sm, sa, se, si);

      r->integrator = REB_INTEGRATOR_MERCURIUS;
      r->dt = 0.1*2.*M_PI;
      r->ri_mercurius.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
      r->ri_mercurius.peri=1;
      r->heartbeat  = heartbeat;

      reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
      e_init = reb_tools_energy(r);
      // system("rm -rf energy.txt");
      // FILE* f = fopen("energy.txt","w");

      reb_integrate(r, 3000 * 2 * M_PI);
      reb_free_simulation(r);
      fprintf(f,"%e,%e\n", peris[i], err);
    }
}


void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}
    if (reb_output_check(r, (0.01 / 365.25) * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        // FILE* f = fopen("energy.txt","a");

        // rotate whole simulation to rotating frame
        //reb_simulation_irotate(r, r1);
        double e = reb_tools_energy(r);
        double err_stage = fabs((e - e_init) / e_init);

        if (err < err_stage){
          err = err_stage;
        }
        // fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e\n",r->t, fabs((e - e_init) / e_init), r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z);
        // fclose(f);

        // reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);
    }
}
