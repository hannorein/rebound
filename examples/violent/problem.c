/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "rebound.h"
#include "integrator_trace.h"
#include "particle.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double tmax = 1e7*2*M_PI;
int nbodies=3;
int first_ejected = 999;
int ind;
int nparticles;

char title[100] = "machine_indep_test_grace_";
char title_stats[100] = "redo_chaotic_bs_time_stats";
char title_remove[100] = "rm -rf machine_indep_test_grace_";

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
    reb_add(r, star);
    double planet_m = 9.55e-4;

    double sma = 5.;
    double add = 0.;
    double delta= 2.5;

    ind = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      strcat(title_remove, argv[1]);
      ind = atoi(argv[1]);
    }

    r->rand_seed = ind;

    for (int i = 0; i < nbodies; i++){
      if (i == 2){
        add = reb_random_uniform(r, -1e-12, 1e-12);
        sma += add;
      }

      reb_add_fmt(r, "m a e inc hash", planet_m, sma, 0.05, (double)i * M_PI / 180., i+1);
      double num = -pow(2., 1./3.) * pow(3., 1./3.) * sma - pow((planet_m / star.m), 1./3.) * delta * sma;
      double denom = -pow(2., 1./3.) * pow(3., 1./3.) + pow((planet_m / star.m), 1./3.) * delta;
      sma = num/denom;
    }

    struct reb_particle* sun = &r->particles[0];
    double min = INFINITY;

    for (int i = 1; i < nbodies+1; i++){
        struct reb_particle* p = reb_get_particle_by_hash(r, i);
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
        if (abs(o.P) < min){
          min = o.P;
        }
    }

    r->integrator = REB_INTEGRATOR_IAS15;
    //r->dt = min * 0.05;//7.108147e-01;//min * 0.010123456;//0.059331635924546614;
    r->ri_ias15.adaptive_mode=2;
    r->track_energy_offset = 1;
    reb_configure_box(r, 10000., 1, 1, 1);
    r->boundary = REB_BOUNDARY_OPEN;
    r->heartbeat  = heartbeat;

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    if (r->heartbeat != NULL){
      system(title_remove);
      FILE* f = fopen(title, "w");
      fprintf(f, "# Seed: %d\n", ind);
      fprintf(f, "t,E,sx,sy,sz");
      for (int i = 1; i < nbodies+1; i++){
        fprintf(f, ",a%d,e%d,i%d,x%d,y%d,z%d",i,i,i,i,i,i);
      }

      fprintf(f, "\n");
      fclose(f);
    }

    e_init = reb_tools_energy(r);
    nparticles = r->N;
    clock_t begin = clock();
    reb_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

/*
    FILE* tf = fopen(title_stats, "a");
    fprintf(tf, "%d,%d,%e,%e,%f,%d\n", ind, surviving_particles, fabs(reb_tools_energy(r) - e_init)/e_init, time_spent, r->t,r->N);
    fclose(tf);
*/
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    // remove particles if needed
    if (r->N < nparticles){
      reb_move_to_com(r);
      nparticles = r->N;
    }

    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, tmax);
    //}

    // Time to first ejection
    // Always track
    /*
    if (first_ejected == 999){ // ejection has not happened yet
      for (unsigned int i = 1; i < nbodies+1; i++){
          struct reb_particle* p = reb_get_particle_by_hash(r, i);
          if (p == NULL){ // first ejected particle
            first_ejected = i;
            printf("First Ejected: %d\n", i);
            break;
          }
      }
      if (first_ejected != 999){ // if detected ejection
        printf("Try to write to file\n");
        FILE* fs = fopen(title_stats, "a");
        printf(title_stats);
        fprintf(fs, "%d,%d,%e\n", ind, first_ejected, r->t);
        fclose(fs);
      }
    }
    */

    if (reb_output_check(r, 10. * 2.*M_PI)){

      FILE* f = fopen(title, "a");
      struct reb_particle* sun = &r->particles[0];
      struct reb_particle com = reb_get_com(r);
      fprintf(f, "%e,%e", r->t, fabs((reb_tools_energy(r) - e_init) / e_init), sun->x,sun->y,sun->z);

      for (unsigned int i = 1; i < nbodies+1; i++){
        struct reb_particle* p = reb_get_particle_by_hash(r, i);
        if (p != NULL){
          struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
          fprintf(f, ",%f,%f,%f,%f,%f,%f", o.a, o.e, o.inc,p->x,p->y,p->z);
        }
        else{
          fprintf(f, ",-1,-1,-1,1e10,1e10,1e10");
        }
      }
      fprintf(f, "\n");
      fclose(f);
    }

}
