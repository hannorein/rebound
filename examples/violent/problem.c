/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
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
double tmax = 1000.*2*M_PI;

char title[100] = "energy.txt";

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
/*
    star.x = -9.184617e-02;
    star.y = 1.183124e+01;
    star.z = 6.603905e-01;
    star.vx = -7.676444e-04;
    star.vy = 2.034279e-04;
    star.vz = 7.427045e-05;
    reb_add(r, star);

    double planet_m = 9.55e-4;
    // double planet_r = 0.000477895;

    struct reb_particle p1 = {0};
    p1.m = planet_m;
    p1.x = -1.066705e+00;
    p1.y = 1.051777e+01;
    p1.z = 4.684022e-01;
    p1.vx = 8.115600e-01;
    p1.vy = -3.859641e-01;
    p1.vz = -5.767662e-02;

    struct reb_particle p2 = {0};
    p2.m = planet_m;
    p2.x = 1.791964e+01;
    p2.y = 1.161478e+01;
    p2.z = 3.701365e+00;
    p2.vx = -8.189499e-03;
    p2.vy = 2.114994e-01;
    p2.vz = -1.787831e-02;

    struct reb_particle p3 = {0};
    p3.m = planet_m;
    p3.x = 7.932107e+01;
    p3.y = -1.241087e+04;
    p3.z = -6.956781e+02;
    p3.vx = 4.455993e-04;
    p3.vy = -3.854883e-02;
    p3.vz = -2.215181e-03;
*/
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
    //reb_add(r, p3);
    reb_add_fmt(r, "m a e inc", planet_m, -742.792028, 5.0, 0.13263699999999998);

    double sma = 5.;
    double delta= 2.5;

    int index = 0;
    if (argc == 2){
      strcat(title, argv[1]);
      index = atoi(argv[1]);
    }
/*
    r->rand_seed = index;

    for (int i = 0; i < nbodies; i++){
      double a_rand = 0.;//reb_random_uniform(r, -1e-10, 1e-10);
      double e_rand = 0.;//reb_random_uniform(r, -1e-10, 1e-10);
      double Om_rand = reb_random_uniform(r, 0, 2 * M_PI);
      double f_rand = reb_random_uniform(r, 0, 2 * M_PI);
      reb_add_fmt(r, "m a e inc Omega f", planet_m, sma, 0.05, (double)i * M_PI / 180., Om_rand, f_rand);
      double num = -pow(2., 1./3.) * pow(3., 1./3.) * sma - pow((planet_m / star.m), 1./3.) * delta * sma;
      double denom = -pow(2., 1./3.) * pow(3., 1./3.) + pow((planet_m / star.m), 1./3.) * delta;

      sma = num/denom;
    }

    struct reb_particle* sun = &r->particles[0];
    double min = INFINITY;

    for (int i = 1; i < nbodies+1; i++){
      struct reb_particle* p = &r->particles[i];
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
      if (abs(o.P) < min){
        min = o.P;
      }
      //printf("%d %f %f\n", i, p->x, p->y);
    }
*/
    r->integrator = REB_INTEGRATOR_TRACE;
    //r->integrator = REB_INTEGRATOR_WHFAST;
    //r->ri_whfast.coordinates = 1;
    r->dt = 7.108147e-01;//min * 0.010123456;//0.059331635924546614;
    r->ri_tr.S_peri = reb_integrator_trace_switch_fdot_peri;
    r->ri_tr.hillfac = 4.;
    //r->ri_tr.peri = 2.;
    r->ri_tr.vfac_p = 16.0;

/*
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->dt = min * 0.010123456;//0.059331635924546614;
    r->ri_mercurius.hillfac = 4.;
*/

    //r->integrator = REB_INTEGRATOR_IAS15;
    r->visualization = REB_VISUALIZATION_NONE;
    r->heartbeat  = heartbeat;

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    if (r->heartbeat != NULL){
      system("rm -rf test_2.txt");
      //FILE* f = fopen(title, "w");
      FILE* f = fopen("test_2.txt", "w");
      fprintf(f, "# Seed: %d\n", index);
      fprintf(f, "t,E,sx,sy,sz");
      for (int i = 1; i < nbodies+1; i++){
        //fprintf(f, ",a%d,x%d,y%d,z%d,vx%d,vy%d,vz%d",i,i,i,i,i,i,i);
        fprintf(f, ",a%d,e%d,i%d,x%d,y%d,z%d",i,i,i,i,i,i);
        //fprintf(f, ",a%d",i);
      }

      fprintf(f, "\n");
      fclose(f);
    }

    tmax =1000;//2e6*2*M_PI; //min * 100000;
    //reb_integrate(r, tmax/2.);
    //r->ri_tr.vfac_p = 1000.0;
    reb_integrate(r, tmax);
/*
    for (int i = 0; i < r->N; i++){
      //fprintf(f, ",a%d,x%d,y%d,z%d,vx%d,vy%d,vz%d",i,i,i,i,i,i,i);
      struct reb_particle* p = &r->particles[i];
      printf("%d %e %e %e %e %e %e\n", i, p->x,p->y,p->z,p->vx,p->vy,p->vz);
      //fprintf(f, ",a%d",i);
    }
    */
    //printf("Total: %d\nInit Peri No Flag: %d\nInit Peri Flag: %d\nNo Flags:%d\nFlagged Peri:%d\nFlagged CE:%d\n", r->ri_tr.delta_Ks[0], r->ri_tr.delta_Ks[1], r->ri_tr.delta_Ks[2], r->ri_tr.delta_Ks[3], r->ri_tr.delta_Ks[4], r->ri_tr.delta_Ks[5]);
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10.*2.*M_PI)){
        reb_output_timing(r, tmax);
    }
    //if (reb_output_check(r, 10. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        //FILE* f = fopen(title, "a");
        FILE* f = fopen("test_2.txt", "a");
        fprintf(f, "%e,%e", r->t, reb_tools_energy(r));
        //fprintf(f, "%f", r->t);
        struct reb_particle* sun = &r->particles[0];

        fprintf(f, ",%e,%e,%e",sun->x,sun->y,sun->z);

        for (int i = 0; i < nbodies; i++){
          struct reb_particle* p = &r->particles[i+1];
          struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
          fprintf(f, ",%f,%f,%f,%f,%f,%f", o.a, o.e, o.inc,p->x,p->y,p->z);
          //printf("%f %f %f\n", o.a, o.e, o.inc);
          //fprintf(f, ",%f,%e,%e,%e,%e,%e,%e", o.a, p->x,p->y,p->z, p->vx, p->vy, p->vz);
          //fprintf(f, ",%e", o.a);
        }

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
        fprintf(f, "\n");
        fclose(f);
    //}
}
