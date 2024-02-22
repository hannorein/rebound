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

double e_start; // initial energy
double tmax = 1e7*2*M_PI;
int nbodies=3;
int first_ejected = 999;
int ind;

char title[100] = "timestamps/ias15_stats_";
char title_stats[100] = "222_ias15_cstats";//"merc_timestamps/mercurius_first_ejection";
char element_stats[100] = "222_ias15_cstats";//"merc_timestamps/mercurius_first_ejection";
char title_remove[100] = "rm -rf timestamps/ias15_stats_";

int main(int argc, char* argv[]){

    // Initialize masses

    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
    reb_simulation_add(r, star);
    double planet_m = 9.55e-4;
    double planet_r = 0.000477895;

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

    add = reb_random_uniform(r, -1e-12, 1e-12);

    double smas[3] = {};

    for (int i = 0; i < nbodies; i++){
      smas[i] = sma;
      reb_simulation_add_fmt(r, "m r a e inc hash", planet_m, planet_r, sma, 0.05, (double)i * M_PI / 180., i+1);
      double num = -pow(2., 1./3.) * pow(3., 1./3.) * sma - pow((planet_m / star.m), 1./3.) * delta * sma;
      double denom = -pow(2., 1./3.) * pow(3., 1./3.) + pow((planet_m / star.m), 1./3.) * delta;

      sma = num/denom;
    }
    // random
    struct reb_particle* pouter = &r->particles[nbodies];
    pouter->x += add;

    struct reb_particle* sun = &r->particles[0];

    double final_a = smas[0] * smas[1] * smas[2] / (smas[0]*smas[1] + smas[1]*smas[2] + smas[0]*smas[2]);
    double final_ts = sqrt(4 * M_PI * M_PI / (r->G * star.m) * final_a * final_a * final_a) / 15.12345;

    reb_simulation_move_to_com(r);

    r->integrator = REB_INTEGRATOR_IAS15;
    //r->ri_trace.peri_crit_distance = 0.25 * final_a;
    //r->ri_trace.S_peri = reb_integrator_trace_switch_peri_distance;
    //r->dt = final_ts;
    r->ri_ias15.adaptive_mode = 2;
    r->exit_max_distance = 1e4;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    //r->ri_trace.peri_mode = REB_TRACE_PERI_FULL_IAS15;
    //r->heartbeat  = heartbeat;               // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    if (r->heartbeat != NULL){
      system(title_remove);
      FILE* f = fopen(title, "w");
      //fprintf(f, "# Seed: %d,%.20e\n", ind, add);
      fprintf(f, "t");
      for (int i = 1; i < nbodies+1; i++){
        fprintf(f, ",a%d,e%d,i%d,x%d,y%d,z%d",i,i,i,i,i,i);
      }

      fprintf(f, "\n");
      fclose(f);
    }


    e_start = reb_simulation_energy(r);
    clock_t begin = clock();

    while (r->t < tmax){
       int retval = reb_simulation_integrate(r, tmax);
       if (retval == REB_STATUS_ESCAPE){
          //printf("Ejecting\n");
          // Track energy offset
          double Ei = reb_simulation_energy(r);

          // Find and remove the particle
          int remove_ind;
          for (int i = 1; i < nbodies+1; i++){

              struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
              //printf("Query %d\n", i);
              if (p != NULL){
                double dx = p->x;
                double dy = p->y;
                double dz = p->z;
                double d2 = dx*dx+dy*dy+dz*dz;

                if (d2>r->exit_max_distance * r->exit_max_distance){
                    remove_ind = i;
                }
              }
          }
          reb_simulation_remove_particle_by_hash(r, remove_ind, 1);
          reb_simulation_move_to_com(r);

          r->energy_offset += Ei - reb_simulation_energy(r);
       }
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    FILE* tf = fopen(title_stats, "a");
    fprintf(tf, "%d,%d,%e,%e\n", ind, r->N-1, fabs((reb_simulation_energy(r) - e_start)/e_start), time_spent);
    fclose(tf);

    struct reb_particle* s = &r->particles[0];
    FILE* ef = fopen(element_stats, "a");
    fprintf(ef, "%d", ind);
    double tot_m = 0.0;
    for (unsigned int i = 1; i < nbodies+1; i++){
        struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
        if (p != NULL){
          struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *s);
          fprintf(ef, ",%f,%f,%f", o.a, o.e, o.inc);
          tot_m += p->m;
        }
        else{
          fprintf(ef, ", , , ");
        }
    }
    int num_collisions = (int) (tot_m / planet_m) - (r->N - 1);
    fprintf(ef, ",%d\n", num_collisions);
    fclose(ef);

    reb_simulation_free(r);
}


void heartbeat(struct reb_simulation* r){
    // remove particles if needed

    //if (reb_simulation_output_check(r, 1000.*2.*M_PI)){
    //    reb_simulation_output_timing(r, tmax);
    //}

    // Time to first ejection
    // Always track
/*
    if (first_ejected == 999){ // ejection has not happened yet
      for (unsigned int i = 1; i < nbodies+1; i++){
          struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
          if (p == NULL){ // first ejected particle
            first_ejected = i;
            //printf("First Ejected: %d\n", i);
            break;
          }
      }
      if (first_ejected != 999){ // if detected ejection
        FILE* fs = fopen(title_stats, "a");
        fprintf(fs, "%d,%d,%e\n", ind, first_ejected, r->t);
        fclose(fs);
        exit(1);
      }
    }
*/
/*
    if (reb_simulation_output_check(r, 10000. * 2.*M_PI)){

      FILE* f = fopen(title, "a");
      struct reb_particle* sun = &r->particles[0];

      //fprintf(f, "%e,%e,%f,%f,%f", r->t, fabs((reb_simulation_energy(r) - e_start) / e_start),sun->x,sun->y,sun->z);
      fprintf(f, "%f", r->t);

      for (unsigned int i = 1; i < nbodies+1; i++){
        struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
        if (p != NULL){
          struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *sun);
          fprintf(f, ",%f,%f,%f,%f,%f,%f", o.a, o.e, o.inc, p->x,p->y,p->z);
        }
        else{
          fprintf(f, ",-1,-1,-1,nan,nan,nan");
        }
      }
      fprintf(f, "\n");
      fclose(f);
    }
*/

}
