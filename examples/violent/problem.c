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

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
    reb_add(r, star);

    r->rand_seed = 1;

    for (int i = 0; i < 10; i++){
      // Add bodies
      double loga = reb_random_uniform(r, log10(0.3), log10(10.));
      reb_add_fmt(r, "m r a f Omega inc", 0.005, 0.0005, pow(10, loga), reb_random_uniform(r, 0., 2.*M_PI), reb_random_uniform(r, 0., 2.*M_PI), (double)i * M_PI / 180.);
    }

    struct reb_particle* sun = &r->particles[0];
    double min = INFINITY;
    for (int i = 1; i < 11; i++){
      struct reb_particle* p = &r->particles[i];
      struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
      if (o.P < min){
        min = o.P;
      }
    }
    r->dt = min * 0.05;//0.059331635924546614;

    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_tr.vSwitch=1;
    r->ri_tr.vfac = 30.;

    //r->ri_tr.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    //r->ri_tr.peri = 0.3;
    r->visualization = REB_VISUALIZATION_NONE;
    r->heartbeat  = heartbeat;

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    system("rm -rf energy.txt");
    FILE* f = fopen("energy.txt","w");
    fprintf(f, "t,E,x0,y0,z0");
    for (int i = 1; i < 11; i++){
      fprintf(f, ",a%d,e%d,i%d,x%d,y%d,z%d",i,i,i,i,i,i);
    }
    fprintf(f, "\n");
    fclose(f);

    reb_integrate(r, 1.*2*M_PI);
    printf("%d %d %d %d %d\n", r->ri_tr.delta_Ks[0], r->ri_tr.delta_Ks[1], r->ri_tr.delta_Ks[2], r->ri_tr.delta_Ks[3], r->ri_tr.delta_Ks[4], r->ri_tr.delta_Ks[5]);
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}
    //if (reb_output_check(r, 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");

        fprintf(f, "%e,%e", r->t, reb_tools_energy(r));

        struct reb_particle* sun = &r->particles[0];
        fprintf(f, ",%e,%e,%e",sun->x,sun->y,sun->z);
        for (int i = 0; i < 10; i++){
          struct reb_particle* p = &r->particles[i+1];
          struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *p, *sun);
          fprintf(f, ",%e,%e,%e,%e,%e,%e", o.a, o.e, o.inc,p->x,p->y,p->z);
        }
        fprintf(f, "\n");
        fclose(f);
    //}
}
