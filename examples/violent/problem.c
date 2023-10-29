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

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;
double emax = 0.0;
double amax = 1.5e4;//1e4;

int nbodies = 3;
double tmax = 1000.*2*M_PI;

char title[100] = "chaotic_mercurius_";
char title_final[100] = "chaotic_mercurius_finals";
char title_remove[100] = "rm -rf chaotic_mercurius_";

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
    reb_add(r, star);
    double planet_m = 9.55e-4;
/*
    star.x = 5.972022e-02;
    star.y = 5.463826e-02;
    star.z = -1.816123e-03;
    star.vx = -2.012882e-04;
    star.vy = 1.092496e-04;
    star.vz = 2.598477e-04;
    reb_add(r, star);

    //double planet_m = 9.55e-4;
    // double planet_r = 0.000477895;

    struct reb_particle p1 = {0};
    p1.m = planet_m;
    p1.x = -1.948851e+01;
    p1.y = -3.835959e+01;
    p1.z = -5.170186e+00;
    p1.vx = 8.342011e-02;
    p1.vy = 1.088821e-01;
    p1.vz = 2.812363e-02;

    struct reb_particle p2 = {0};
    p2.m = planet_m;
    p2.x = -9.339950e-01;
    p2.y = 1.092101e+00;
    p2.z = 3.760712e+00;
    p2.vx = 9.243832e-02;
    p2.vy = -8.184879e-02;
    p2.vz = -2.935659e-01;

    struct reb_particle p3 = {0};
    p3.m = planet_m;
    p3.x = -4.211175e+01;
    p3.y = -1.994535e+01;
    p3.z = 3.311174e+00;
    p3.vx = 3.491459e-02;
    p3.vy = -1.414308e-01;
    p3.vz = -6.649501e-03;

    // In DH!
    star.x = 0.0;
    star.y = 0.0;
    star.z = 0.0;
    star.vx = -7.676444e-04;
    star.vy = 2.034279e-04;
    star.vz = 7.427045e-05;
    reb_add(r, star);

    double planet_m = 9.55e-4;
    // double planet_r = 0.000477895;

    struct reb_particle p1 = {0};
    p1.m = planet_m;
    p1.x = -9.748588e-01;
    p1.y = -1.313470e+00;
    p1.z = -1.919883e-01;
    p1.vx = 8.115600e-01;
    p1.vy = -3.859641e-01;
    p1.vz = -5.767662e-02;

    struct reb_particle p2 = {0};
    p2.m = planet_m;
    p2.x = 1.801149e+01;
    p2.y = -2.164600e-01;
    p2.z = 3.040974e+00;
    p2.vx = -8.189499e-03;
    p2.vy = 2.114994e-01;
    p2.vz = -1.787831e-02;

    struct reb_particle p3 = {0};
    p3.m = planet_m;
    p3.x = 7.941292e+01;
    p3.y = -1.242270e+04;
    p3.z = -6.963385e+02;
    p3.vx = 4.455993e-04;
    p3.vy = -3.854883e-02;
    p3.vz = -2.215181e-03;


    reb_add(r, p1);
    reb_add(r, p2);
    reb_add(r, p3);
    //reb_add_fmt(r, "m a e inc", planet_m, -742.792028, 5.0, 0.13263699999999998);
*/
    double sma = 5.;
    double delta= 2.5;

    int index = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      strcat(title_remove, argv[1]);
      index = atoi(argv[1]);
    }

    r->rand_seed = index;

    for (int i = 0; i < nbodies; i++){
      double Om_rand = 0.0;//reb_random_uniform(r, 0, 2 * M_PI);
      double f_rand = 0.0;//reb_random_uniform(r, 0, 2 * M_PI);
      if (i == 2){
        sma += reb_random_uniform(r, -1e-12, 1e-12);
      }
      reb_add_fmt(r, "m a e inc Omega f", planet_m, sma, 0.05, (double)i * M_PI / 180., Om_rand, f_rand);
      double num = -pow(2., 1./3.) * pow(3., 1./3.) * sma - pow((planet_m / star.m), 1./3.) * delta * sma;
      double denom = -pow(2., 1./3.) * pow(3., 1./3.) + pow((planet_m / star.m), 1./3.) * delta;

      sma = num/denom;
    }

    struct reb_particle* sun = &r->particles[0];
    double min = INFINITY;

    for (int i = 1; i < nbodies+1; i++){
      struct reb_particle* p = &r->particles[i];
      p->hash = i;
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
      if (abs(o.P) < min){
        min = o.P;
      }
      //printf("%d %f %f\n", i, p->x, p->y);
    }

    //struct reb_particle* sun = &r->particles[0];
    //r->integrator = REB_INTEGRATOR_WHFAST;
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    //r->ri_ias15.adaptive_mode = 2;
    //r->ri_whfast.coordinates = 1;
    r->dt = min * 0.05; //7.108147e-01;//min * 0.010123456;//0.059331635924546614;
    //r->ri_tr.S_peri = reb_integrator_trace_switch_fdot_peri;
    r->ri_mercurius.hillfac = 4.;
    //r->ri_tr.peri = 2.;
    //r->ri_tr.vfac_p = 32.0;

/*
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->dt = min * 0.010123456;//0.059331635924546614;
    r->ri_mercurius.hillfac = 4.;
*/

    //r->integrator = REB_INTEGRATOR_BS;
    r->visualization = REB_VISUALIZATION_NONE;
    r->heartbeat  = heartbeat;

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    if (r->heartbeat != NULL){
      system(title_remove);
      FILE* f = fopen(title, "w");
      //FILE* f = fopen("test.txt", "w");
      //fprintf(f, "# Seed: %d\n", index);
      fprintf(f, "t,E");
      for (int i = 1; i < nbodies+1; i++){
        //fprintf(f, ",a%d,x%d,y%d,z%d,vx%d,vy%d,vz%d",i,i,i,i,i,i,i);
        //fprintf(f, ",hash%d,a%d,e%d,i%d,x%d,y%d,z%d",i,i,i,i,i,i,i);
        fprintf(f, ",hash%d,d%d,ej%d",i,i,i);
      }

      fprintf(f, "\n");
      fclose(f);
    }

    tmax = 1e7 * 2 * M_PI;//1e6*2*M_PI;//2e6*2*M_PI; //min * 100000;
    e_init = reb_tools_energy(r);
    //clock_t begin = clock();
    reb_integrate(r, tmax);

    FILE* tf = fopen(title_final, "a");
    fprintf(tf, "%d", index);
    for (int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
      fprintf(tf, ",%u,%e,%e,%e", p->hash, o.a, o.e, o.inc);
      //fprintf(f, ",%e", o.a);
    }
    fprintf(tf, "\n");
    fclose(tf);
    // clock_t end = clock();

    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //FILE* f = fopen(title, "a");
    //fprintf(f, "%f\n", time_spent);
    //fclose(f);
    //printf("%f\n", time_spent);

    //for (int i = 0; i < r->N; i++){
      //fprintf(f, ",a%d,x%d,y%d,z%d,vx%d,vy%d,vz%d",i,i,i,i,i,i,i);
    //  struct reb_particle* p = &r->particles[i];
    //  printf("%d %e %e %e %e %e %e\n", i, p->x,p->y,p->z,p->vx,p->vy,p->vz);
      //fprintf(f, ",a%d",i);
    //}

    //printf("Total: %d\nInit Peri No Flag: %d\nInit Peri Flag: %d\nNo Flags:%d\nFlagged Peri:%d\nFlagged CE:%d\n", r->ri_tr.delta_Ks[0], r->ri_tr.delta_Ks[1], r->ri_tr.delta_Ks[2], r->ri_tr.delta_Ks[3], r->ri_tr.delta_Ks[4], r->ri_tr.delta_Ks[5]);
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    // remove particles if needed
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, tmax);
    //}
    if (reb_output_check(r, 1. * 2.*M_PI)){
      FILE* f = fopen(title, "a");
      fprintf(f, "%e,%e", r->t, fabs((reb_tools_energy(r) - e_init) / e_init));
      struct reb_particle* sun = &r->particles[0];
      for (unsigned int i = 1; i < r->N; i++){
        unsigned int ej = 0;
        struct reb_particle* p = &r->particles[i];
        double dx = p->x - sun->x;
        double dy = p->y - sun->y;
        double dz = p->z - sun->z;
        double d = sqrt(dx*dx + dy*dy + dz*dz);
        //struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
        fprintf(f, ",%u,%f", p->hash, d);
        if (d > amax){
          //printf("Ejection %u\n", i);
          reb_remove(r, i, 1);
          e_init = reb_tools_energy(r);
          ej = 1;
        }
        //printf("Here\n");
        fprintf(f, ",%u", ej);
      }
      fprintf(f, "\n");
      fclose(f);
    }
    //if (reb_output_check(r, 1. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        //FILE* f = fopen(title, "a");
        //fprintf(f, "%f", r->t);
        //struct reb_particle* sun = &r->particles[0];

        // fprintf(f, ",%e,%e,%e",sun->x,sun->y,sun->z);
/*
        for (int i = 1; i < r->N; i++){
          struct reb_particle* p = &r->particles[i];
          struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
          //fprintf(f, ",%u,%f,%f,%f,%f,%f,%f", p->hash, o.a, o.e, o.inc,p->x,p->y,p->z);
          //printf("%f %f %f\n", o.a, o.e, o.inc);
          fprintf(f, ",%f,%e,%e,%e,%e,%e,%e", o.a, p->x,p->y,p->z, p->vx, p->vy, p->vz);
          //fprintf(f, ",%e", o.a);
        }
*/

/*
        struct reb_particle* p1 = &r->particles[1];
        struct reb_particle* p2 = &r->particles[2];
        const double dx = p1->x - p2->x;
        const double dy = p1->y - p2->y;
        const double dz = p1->z - p2->z;
        const double dij = sqrt(dx*dx + dy*dy + dz*dz);

        const double dvx = p1->vx - p2->vx;
        const double dvy = p1->vy - p2->vy;
        const double dvz = p1->vz - p2->vz;
        const double dvij = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

        fprintf(f, "%e,%e,%e,%e", r->t, reb_tools_energy(r), dij, dvij);
*/
        //fprintf(f, "\n");
        //fclose(f);
  //  }
}