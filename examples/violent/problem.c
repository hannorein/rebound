/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "integrator_trace.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;
double emax = 0.0;

// Values from Issue 544
double as[10] = {0.4805534574239261,
  0.7337708461939174,
  2.9472592774364457,
  0.33135308945198694,
  4.344916055712149,
  3.7665109455686245,
  7.077599283999203,
  2.002750283874912,
  0.6411739552952435,
  0.6527258031474394};

double fs[10] = {3.141592653589793, 2.940157490939905, 1.5707963267948966, 1.5707963267948966, 2.9083687607633983,
  6.23982548472831, 0.15988105992264678, 0.12859413226990846, 6.2831853071795845, 6.104371966391728};

double Omegas[10] = {0.0, 3.1129104598773245, -1.3274903358954842, -1.0319182859768408, 0.01323272347183684, 1.4373554275242735,
  0.1922025319432963, -0.38233708635963903, 2.6522366656182905, 2.751328894621499};

int nbodies = 3;
double tmax = 200000.*2*M_PI;

char title[100] = "violent_ias15_";

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
    reb_add(r, star);

    double planet_m = 9.55e-4;
    double planet_r = 0.000477895;

    double sma = 1.;
    double delta= 3.;

    int index = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      index = atoi(argv[1]);
    }

    r->rand_seed = index;

    for (int i = 0; i < nbodies; i++){
      double a_rand = reb_random_uniform(r, -1e-10, 1e-10);
      double e_rand = reb_random_uniform(r, -1e-10, 1e-10);
      double Om_rand = reb_random_uniform(r, 0, 1e-10);
      double f_rand = reb_random_uniform(r, 0, 1e-10);
      reb_add_fmt(r, "m r a e inc Omega f", planet_m, planet_r, sma + a_rand, 0.01 + e_rand, (double)i * M_PI / 180., Om_rand, f_rand);
      double num = -pow(2., 1./3.) * pow(3., 1./3.) * sma - pow((planet_m / star.m), 1./3.) * delta * sma;
      double denom = -pow(2., 1./3.) * pow(3., 1./3.) + pow((planet_m / star.m), 1./3.) * delta;

      sma = num/denom;
    }

    struct reb_particle* sun = &r->particles[0];
    double min = INFINITY;
    for (int i = 1; i < nbodies+1; i++){
      struct reb_particle* p = &r->particles[i];
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
      if (o.P < min){
        min = o.P;
      }
    }

/*
    r->integrator = REB_INTEGRATOR_TRACE;
    r->dt = min * 0.05;//0.059331635924546614;
    r->ri_tr.S = reb_integrator_trace_switch_velocity;
    r->ri_tr.vfac = 5.;
    r->ri_tr.vfac_p = 16.;
    r->ri_tr.ats = 0;
*/


    r->integrator = REB_INTEGRATOR_IAS15;
    r->visualization = REB_VISUALIZATION_NONE;
    r->heartbeat  = heartbeat;

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    if (r->heartbeat != NULL){
      //system("rm -rf energy.txt");
      FILE* f = fopen(title, "w");
      fprintf(f, "# Seed: %d\n", index);
      fprintf(f, "t,E");
      for (int i = 1; i < nbodies+1; i++){
        // fprintf(f, ",a%d,e%d,i%d,x%d,y%d,z%d",i,i,i,i,i,i);
        fprintf(f, ",a%d",i);
      }
      fprintf(f, "\n");
      fclose(f);
    }


    reb_integrate(r, tmax);
    //reb_step(r);
    //printf("Total: %d\nInit Peri No Flag: %d\nInit Peri Flag: %d\nNo Flags:%d\nFlagged Peri:%d\nFlagged CE:%d\n", r->ri_tr.delta_Ks[0], r->ri_tr.delta_Ks[1], r->ri_tr.delta_Ks[2], r->ri_tr.delta_Ks[3], r->ri_tr.delta_Ks[4], r->ri_tr.delta_Ks[5]);
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 1.*2.*M_PI)){
    //    reb_output_timing(r, tmax);
    //}
    //if (reb_output_check(r, 100. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        FILE* f = fopen(title, "a");

        fprintf(f, "%e,%e", r->t, reb_tools_energy(r));

        struct reb_particle* sun = &r->particles[0];
        //fprintf(f, ",%e,%e,%e",sun->x,sun->y,sun->z);
        for (int i = 0; i < nbodies; i++){
          struct reb_particle* p = &r->particles[i+1];
          struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
          //fprintf(f, ",%e,%e,%e,%e,%e,%e", o.a, o.e, o.inc,p->x,p->y,p->z);
          fprintf(f, ",%e", o.a);
        }
        fprintf(f, "\n");
        fclose(f);
    //}
}
