/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e1_init; // initial energy
double e2_init;
double etot = 0.0;
double emax = 0.0;

double jacobi_dh(struct reb_simulation* r, int p){
  // Jacobi constant for the restricted 3-body problem
  // r needs to be interpreted in the rotating frame
  struct reb_particle* sun = &r->particles[0];
  struct reb_particle* jup = &r->particles[1];
  struct reb_particle* test = &r->particles[p];

  struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *jup, *sun);

  struct reb_vec3d x = {.x = test->x, .y=test->y, .z=test->z};
  struct reb_vec3d v = {.x = test->vx, .y=test->vy, .z=test->vz};

  double kinetic = 0.5 * reb_vec3d_length_squared(v);

  double d0 = sqrt(reb_vec3d_length_squared((struct reb_vec3d){.x = test->x - sun->x, .y = test->y - sun->y, .z =test->z - sun->z}));
  double d1 = sqrt(reb_vec3d_length_squared((struct reb_vec3d){.x = test->x - jup->x, .y = test->y - jup->y, .z =test->z - jup->z}));

  double U = -((r->G * sun->m) / d0) - ((r->G * jup->m) / d1);

  double omega = sqrt((r->G * (sun->m + jup->m) / (o.a * o.a * o.a)));
  double L =  x.x * v.y - x.y * v.x;

  double jac = kinetic + U - omega * L;
  return jac;

}


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    //r->dt = (8./365.)*2.*M_PI;

    // Command line arguments for new timestep
    //double new_dt = 0.0;
    //if (argc == 2){
    //  new_dt = (atof(argv[1]) * 1.0);
    r->dt = (8. / 365.) * 2. * M_PI;
    //}
    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_trace.hillfac = 4.;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    r->ri_trace.peri = 0.01;
    r->N_active = 2;
    r->visualization= REB_VISUALIZATION_NONE;
    //r->heartbeat  = heartbeat;
    //r->ri_mercurius.hillfac = 4;

    // Initialize masses
    struct reb_particle star = {0};
    star.m = 1;
    struct reb_particle jup = {0};
    jup.m = 0.01 / (star.m - 0.01);

    // velocities
    double a = 5.2;
    double e = 0.0;
    star.x = -(jup.m / (star.m + jup.m)) * (a * (1 + e));
    star.vy = -(jup.m / (star.m + jup.m)) * sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
    reb_add(r, star);

    jup.x = (star.m / (star.m + jup.m)) * (a * (1 + e));
    jup.vy = (star.m / (star.m + jup.m)) * sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
    reb_add(r, jup);

    // Test particle
    struct reb_particle test1 = {0};
    double xhel1 = 4.42;
    double vhel1 = 0.0072 * (365.25) * (1 / (2 * M_PI)); // days to REBOUND years
    test1.x = xhel1 + star.x;
    test1.vy = vhel1 + star.vy;

    // 2nd test particle slightly displaced
    struct reb_particle test2 = {0};
    double xhel2 = 4.41;
    double vhel2 = 0.0072 * (365.25) * (1 / (2 * M_PI)); // days to REBOUND years
    test2.x = xhel2 + star.x;
    test2.vy = vhel2 + star.vy;

    reb_add(r, test1);
    reb_add(r, test2);

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e1_init = jacobi_dh(r, 2);
    e2_init = jacobi_dh(r, 3);
    //system("rm -rf energy.txt");

    reb_integrate(r, 500000.*11.86*2.*M_PI);
    //reb_integrate(r, 120.);
    //reb_integrate(r, 119.);
    reb_free_simulation(r);
    // printf("%f %e %e\n", new_dt, etot, emax);
}

void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}
    double e1 = jacobi_dh(r, 2);
    double e2 = jacobi_dh(r, 3);

    double estep = fabs((e1-e1_init)/e1_init + (e2-e2_init)/e2_init);
    if (estep > emax){
      emax = estep;
    }

    etot += ((e1-e1_init)/e1_init + (e2-e2_init)/e2_init);
    //if (reb_output_check(r, (40. / 365.25) * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
    //    FILE* f = fopen("energy.txt","a");

        // rotate whole simulation to rotating frame
        //struct reb_vec3d v1 = {.x = r->particles[1].x, .y = r->particles[1].y, .z = r->particles[1].z};
        //struct reb_vec3d v2 = {.x = (1 / (1 + (0.01 / (1 - 0.01)))) * 5.2, .y = 0, .z = 0};
        //struct reb_rotation r1 = reb_rotation_init_from_to(v1, v2);
        //reb_simulation_irotate(r, r1);
        //reb_integrator_synchronize(r);
    //    fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",r->t, (e1-e1_init)/e1_init, (e2-e2_init)/e2_init, r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z, r->particles[3].x, r->particles[3].y, r->particles[3].z);
    //    fclose(f);

        //reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);

        //struct reb_rotation inverse = reb_rotation_inverse(r1);
        //reb_simulation_irotate(r, inverse);
    //}
}
