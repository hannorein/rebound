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
int Nparticles = 10000;
double tmax = 1e3;

char title[100] = "ias15_accretion";
char title_remove[100] = "rm -rf ias15_accretion";

char remove_snapshots[100] = "rm -rf *kokubo_trace_snapshot_*";
char snapshot_1[100] = "kokubo_trace_snapshot_1";
char snapshot_2[100] = "kokubo_trace_snapshot_2";
char snapshot_3[100] = "kokubo_trace_snapshot_3";
char snapshot_4[100] = "kokubo_trace_final_orbit_snapshot_4";

double snap1_time = 0.0;
double snap2_time = 30.0 * 2 * M_PI;
double snap3_time = 100.0 * 2 * M_PI;
double snap4_time = 1000.0 * 2 * M_PI;
int snap1=1;
int snap2=1;
int snap3=1;
// accretion of the moon
int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    //r->ri_ias15.adaptive_mode=2;

    r->dt = 0.1;//0.003 / 10.56789;
    //r->integrator = REB_INTEGRATOR_IAS15;
    //r->softening = 3e-8;
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = reb_collision_resolve_merge;

    // Initialize masses
    struct reb_particle earth = {0};
    earth.m = 1;
    double earth_r = 1./2.9;//0.00592; // in units of Roche radius
    earth.r = earth_r;
    reb_simulation_add(r, earth);

    r->rand_seed = 1;
    r->heartbeat = heartbeat;

    double lunar_mass = 0.0123;

    double mtot = 0.;
    double rad_tot = 0;
    // Test particles
    for (unsigned int i = 0; i < Nparticles; i++){
      double m = reb_random_powerlaw(r, 1e-5 * lunar_mass, 1e-5 * lunar_mass * 1000., -1.5);
      mtot += m;
      double rad = pow(m, 1./3.) * 1.185 * earth_r;
      rad_tot += rad;
      double a = reb_random_powerlaw(r, earth_r, 1.0, -1.);
      double e = reb_random_rayleigh(r, 0.3);
      if (e > 1.){
        e = 0.99;
      }
      double inc = reb_random_rayleigh(r, 0.6);
      double Omega = reb_random_uniform(r, 0, 2 * M_PI);
      double omega = reb_random_uniform(r, 0, 2 * M_PI);
      double f = reb_random_uniform(r, 0, 2 * M_PI);
      reb_simulation_add_fmt(r, "primary m r a e inc Omega omega f", earth, m, rad, a, e, inc, Omega, omega, f);
    }
    printf("%f\n", mtot/lunar_mass);
    printf("%f\n", rad_tot/Nparticles);
    //exit(1);
    system(title_remove);
    system(remove_snapshots);
    //exit(1);

    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.


    FILE* f = fopen(snapshot_1,"a");
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];

      double dx = p->x - earth.x;
      double dy = p->y - earth.y;
      double dz = p->z - earth.z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f, "%f,%f,%f\n",p->m,r,dz);
    }
    fclose(f);


    e_init = reb_simulation_energy(r);
    //printf("Initial Energy: %e\n", e_init);
    clock_t begin = clock();
    reb_simulation_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("Final Energy: %e\n", reb_simulation_energy(r));
    //printf("Conservation: %f\n", (reb_simulation_energy(r) - e_init) / e_init);
    //printf("Final N: %d\n", r->N);
    printf("Time Spent: %f\n", time_spent);

    FILE* f4 = fopen(snapshot_4,"a");
    struct reb_particle* e = &r->particles[0];
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];
      struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *e);

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f4, "%f,%f,%f,%f,%f,%f\n",p->m,r,dz,o.a,o.e,o.inc);
    }
    fclose(f4);

    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
  if (reb_simulation_output_check(r, 0.01*2.*M_PI)){
      reb_simulation_output_timing(r, tmax);
  }

  if (snap2 && r->t > snap2_time){
    struct reb_particle* e = &r->particles[0];
    snap2 = 0;
    FILE* f = fopen(snapshot_2,"a");
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f, "%f,%f,%f\n",p->m,r,dz);
    }
    fclose(f);
  }

  if (snap3 && r->t > snap3_time){
    struct reb_particle* e = &r->particles[0];
    snap3 = 0;
    FILE* f = fopen(snapshot_3,"a");
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f, "%f,%f,%f\n",p->m,r,dz);
    }
    fclose(f);
  }

  //if (r->t < 10. * 2 * M_PI){

  if (reb_simulation_output_check(r, 0.1)){
    FILE* f = fopen(title,"a");
    //fprintf(f, "%f,%e,%d\n", r->t, fabs((reb_tools_energy(r) - e_init) / e_init), r->N-1);
    fprintf(f, "%f,%d\n", r->t, r->N-1);
    fclose(f);
  }

  //}
  /*
  if (r->t > 10. * 2 * M_PI && r->t < 10. * 2 * M_PI){
    if (reb_output_check(r, 10. * 2.*M_PI)){
      FILE* f = fopen(title,"a");
      fprintf(f, "%f,%e,%d\n", r->t, fabs((reb_tools_energy(r) - e_init) / e_init), r->N-1);
      fclose(f);
    }
  }

  if (r->t > 100. * 2 * M_PI){
    if (reb_output_check(r, 100. * 2.*M_PI)){
      FILE* f = fopen(title,"a");
      fprintf(f, "%f,%e,%d\n", r->t, fabs((reb_tools_energy(r) - e_init) / e_init), r->N-1);
      fclose(f);
    }
  }
  */
}
