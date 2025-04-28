#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"

char TITLE[100] = "IAS15_merge_simarchive_";
char STATS[100] = "IAS15_merge_stats";

int trace_close_encounters = 0;
int mercurius_close_encounters = 0;

double RCRIT = 100000.;
int IND;


#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    // Returns the maximum of a and b

int constant_trace_switch(struct reb_simulation* const r, const unsigned int i, const unsigned int j){
  const double dxi  = r->particles[i].x;
  const double dyi  = r->particles[i].y;
  const double dzi  = r->particles[i].z;

  const double dxj  = r->particles[j].x;
  const double dyj  = r->particles[j].y;
  const double dzj  = r->particles[j].z;

  const double dx = dxi - dxj;
  const double dy = dyi - dyj;
  const double dz = dzi - dzj;
  const double rp = dx*dx + dy*dy + dz*dz;

  double dcriti = 0.;
  double dcritj = 0.;
  if (i == 0){
    dcriti = r->particles[i].r;
  }
  else{
    dcriti = r->particles[i].m * RCRIT;
    dcritj = r->particles[j].m * RCRIT;
  }

  double dcritmax = MAX(dcriti, dcritj);

  return rp<dcritmax*dcritmax;
}

int lol(struct reb_simulation* const r, const unsigned int j){
  return 1;
}

double get_radii(double m, double rho){
    return pow((3.*m)/(4.*M_PI*rho),1./3.);
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 39.476926421373;
    //r->exit_max_distance = 100.;

    // The random seed is passed as a command line argument
    if (argc == 2){
      IND = atoi(argv[1]);
      strcat(TITLE, argv[1]);
    }
    r->rand_seed = IND;

    // MERCURIUS settings
    // ---------------------------------------------------------
    //r->integrator = REB_INTEGRATOR_MERCURIUS;
    //r->dt = 5./365.;
    // ---------------------------------------------------------

    // TRACE settings
    // ---------------------------------------------------------
    //r->integrator = REB_INTEGRATOR_TRACE;
    //r->ri_trace.S = constant_trace_switch;
    //r->ri_trace.S_peri = reb_integrator_trace_switch_peri_none;
    //r->dt = 5./365.;
    // --------------------------------------------------------

    // IAS15 settings
    // ---------------------------------------------------------
    r->integrator = REB_INTEGRATOR_IAS15;
    // --------------------------------------------------------
    r->collision = REB_COLLISION_DIRECT;
    //r->collision_resolve = reb_collision_resolve_fragment;
    r->collision_resolve = reb_collision_resolve_merge;

    //Assigning mass and number of planetary embryos and planetesimals
    struct reb_particle star = {0};
    star.m = 1.0;
    star.r = 0.00465;
    reb_simulation_add(r, star);

    // Constants for mass range
    double lunar_mass = 3.8e-8;
    double earth_mass = 3e-6;
    double mass_min = 0.6 * lunar_mass;
    double mass_max = 0.2 * earth_mass;
    double rho = 5.05e6; //3 g/cm^3

    // Add 30 planetary embryos
    for (int i = 0; i < 30; i++) {
        double a = reb_random_uniform(r, 0.1, 0.5);     // semi-major axis in AU
        double e = reb_random_uniform(r, 0.0, 0.01);    // eccentricity
        double inc = reb_random_uniform(r, 0.0, M_PI/180.);                        // inclination
        double omega = reb_random_uniform(r, 0.0, 2.*M_PI); // argument of periapsis
        double Omega = reb_random_uniform(r, 0.0, 2.*M_PI); // longitude of ascending node
        double f = reb_random_uniform(r, 0.0, 2.*M_PI);     // mean anomaly
        double m = reb_random_uniform(r, mass_min, mass_max);  // in solar masses

        struct reb_particle emb = reb_particle_from_orbit(r->G, star, m, a, e, inc, Omega, omega, f);
        emb.r = get_radii(m, rho);

        reb_simulation_add(r, emb);
    }


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 1e2;
    //reb_simulation_save_to_file_interval(r,TITLE,1.e2);
    reb_simulation_integrate(r, run_time);

    FILE* stats = fopen(STATS,"a");
    fprintf(stats,"Seed:%d Total Steps: %d Rejected Steps: %ld MERCRUIUS Encounters: %d TRACE encounters: %d\n", IND, r->steps_done, r->ri_trace.step_rejections, mercurius_close_encounters, trace_close_encounters);
    fclose(stats);

    reb_simulation_free(r);
}

//void heartbeat(struct reb_simulation* r){
//    if (reb_simulation_output_check(r, 1.e3)){
//        reb_simulation_output_timing(r, 0);
//        printf("Walltime(s) = %f \n", r->walltime);
//    }
//}
