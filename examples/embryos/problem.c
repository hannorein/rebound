#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rebound.h"

double min_frag_mass = 1.4e-8;

char TITLE[100] = "simarchive_fragment2_";
char STATS[100] = "stats_fragment2_";

int trace_close_encounters = 0;
int mercurius_close_encounters = 0;
double RHO = 5.05e6; //3 g/cm^3

double RCRIT = 100000.;
int IND;

int cols = 0;


#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    // Returns the maximum of a and b

double get_radii(double m, double rho){
    return pow((3.*m)/(4.*M_PI*rho),1./3.);
}

int fragment(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].last_collision==r->t || r->particles[c.p2].last_collision==r->t) return 0;

    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    unsigned int i = c.p1;
    unsigned int j = c.p2;   //want j to be removed particle
    if (j<i){
        swap = 1;
        i = c.p2;
        j = c.p1;
    }

    if (i == 0){
      // Collision with central body, just remove
      return swap?1:2;
    }

    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
    struct reb_particle com = reb_particle_com_of_pair(*pi, *pj);
    double total_mass = pi->m + pj->m;
    double total_rad = pi->r + pj->r;
    double rfactor = 10; // avoid collision cascades

    // overwrite particle 1 with large fragment
    double vesc = sqrt(2.*r->G*total_mass);
    struct reb_vec3d vcom = {.x=com.vx,.y=com.vy,.z=com.vz};
    struct reb_vec3d vcom_hat = reb_vec3d_normalize(vcom);

    pi->m = total_mass/2.;
    pi->r = get_radii(pi->m, RHO);
    pi->x = com.x + 4.*rfactor*total_rad*vcom_hat.x;
    pi->y = com.y + 4.*rfactor*total_rad*vcom_hat.y;
    pi->z = com.z + 4.*rfactor*total_rad*vcom_hat.z;
    pi->vx = com.vx + 2.*vesc*vcom_hat.x;
    pi->vy = com.vy + 2.*vesc*vcom_hat.y;
    pi->vz = com.vz + 2.*vesc*vcom_hat.z;

    if (total_mass/4. > min_frag_mass){
      // Two fragments in the opposite direction
      struct reb_particle frag1 = {0};
      frag1.m = total_mass/4.;
      frag1.r = get_radii(frag1.m,RHO);
      frag1.x = com.x - 6.*rfactor*total_rad*vcom_hat.x;
      frag1.y = com.y - 6.*rfactor*total_rad*vcom_hat.y;
      frag1.z = com.z - 6.*rfactor*total_rad*vcom_hat.z;
      frag1.vx = com.vx - 2.*vesc * vcom_hat.x;
      frag1.vy = com.vy - 2.*vesc * vcom_hat.y;
      frag1.vz = com.vz - 2.*vesc * vcom_hat.z;

      struct reb_particle frag2 = {0};
      frag2.m = total_mass/4.;
      frag2.r = get_radii(frag2.m,RHO);
      frag2.x = com.x - 2.*rfactor*total_rad*vcom_hat.x;
      frag2.y = com.y - 2.*rfactor*total_rad*vcom_hat.y;
      frag2.z = com.z - 2.*rfactor*total_rad*vcom_hat.z;
      frag2.vx = com.vx - 2.*vesc*vcom_hat.x;
      frag2.vy = com.vy - 2.*vesc*vcom_hat.y;
      frag2.vz = com.vz - 2.*vesc*vcom_hat.z;

      reb_simulation_add(r,frag1);
      reb_simulation_add(r,frag2);
    }
    return swap?1:2;
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.)){
        //reb_simulation_output_timing(r, 0);
        printf("Walltime(s) = %f %d\n", r->walltime, r->N);
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 39.476926421373;
    //r->exit_max_distance = 100.;

    // The random seed is passed as a command line argument
    int integrator = 0;
    if (argc == 3){
      integrator = atoi(argv[1]);
      IND = atoi(argv[2]);
    }
    r->rand_seed = IND;

    switch(integrator){
      case 0:
      {
        printf("Didn't Set Integrator!");
        exit(1);
      }
      break;
      case 1:
      {
        //printf("IAS15\n");
        r->integrator = REB_INTEGRATOR_IAS15;
        strcat(TITLE, "IAS15_");
        strcat(STATS, "IAS15");
      }
      break;
      case 2:
      {
        //printf("BS\n");
        r->integrator = REB_INTEGRATOR_BS;
        strcat(TITLE, "BS_");
        strcat(STATS, "BS");
      }
      break;
      case 3:
      {
        r->integrator = REB_INTEGRATOR_MERCURIUS;
        r->dt = 5./365.;
        //printf("MERCURIUS\n");
        strcat(TITLE, "MERCURIUS_");
        strcat(STATS, "MERCURIUS");
      }
      break;
      case 4:
      {
        r->integrator = REB_INTEGRATOR_TRACE;
        //r->ri_trace.S = constant_trace_switch;
        r->ri_trace.S_peri = reb_integrator_trace_switch_peri_none;
        r->ri_trace.r_crit_hill = 3.*1.21;
        r->dt = 5./365.;
        //printf("TRACE\n");
        strcat(TITLE, "TRACE_");
        strcat(STATS, "TRACE");
      }
      break;
    }
    strcat(TITLE, argv[2]);
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = fragment;
    //r->heartbeat = heartbeat;

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
        emb.r = get_radii(m, RHO);

        reb_simulation_add(r, emb);
    }


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    double run_time = 1e5;
    //printf("%d\n",r->N);
    reb_simulation_save_to_file_interval(r,TITLE,1.e2);
    reb_simulation_integrate(r, run_time);
    //printf("%d\n",r->N);

    FILE* stats = fopen(STATS,"a");
    fprintf(stats,"Seed:%d Total Steps: %d Rejected Steps: %ld MERCRUIUS Encounters: %d TRACE encounters: %d\n", IND, r->steps_done, r->ri_trace.step_rejections, mercurius_close_encounters, trace_close_encounters);
    fclose(stats);

    reb_simulation_free(r);
}
