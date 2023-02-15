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
/*
double jacobi(struct reb_simulation* r){
  // Jacobi constant for the restricted 3-body problem
  // r needs to be interpreted in the rotating frame
  struct reb_particle* sun = &r->particles[0];
  struct reb_particle* jup = &r->particles[1];
  struct reb_particle* test = &r->particles[2];

  struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *jup, *sun);

  struct reb_vec3d x = {.x = test->x, .y=test->y, .z=test->z};
  struct reb_vec3d v = {.x = test->vx, .y=test->vy, .z=test->vz};

  double kinetic = 0.5 * reb_vec3d_length_squared(v);

  double d0 = sqrt(reb_vec3d_length_squared((struct reb_vec3d){.x = test->x - sun->x, .y = test->y - sun->y, .z =test->z - sun->z}));
  double d1 = sqrt(reb_vec3d_length_squared((struct reb_vec3d){.x = test->x - jup->x, .y = test->y - jup->y, .z =test->z - jup->z}));

  struct reb_vec3d omega = {.x=0, .y=0, .z= (r->G * (sun->m + jup->m) / (o.a * o.a * o.a))};

  double phi_eff = -((r->G * sun->m) / d0) - ((r->G * jup->m) / d1) - 0.5 * reb_vec3d_length_squared(reb_vec3d_cross(omega, x));
  return kinetic + phi_eff;

}
*/
double jacobi_dh(struct reb_simulation* r){
  // Jacobi constant for the restricted 3-body problem
  // r needs to be interpreted in the rotating frame
  struct reb_particle* sun = &r->particles[0];
  struct reb_particle* jup = &r->particles[1];
  struct reb_particle* test = &r->particles[2];

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
    r->dt = 0.1*2.*M_PI;
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->ri_mercurius.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    r->heartbeat  = heartbeat;

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
    struct reb_particle test = {0};
    double xhel = 4.42;
    double vhel = 0.0072 * (365.25) * (1 / (2 * M_PI)); // days to REBOUND years

    test.x = xhel + star.x;
    test.vy = vhel + star.vy;
    reb_add(r, test);


    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_tools_energy(r);
    system("rm -rf energy.txt");

    reb_integrate(r, 50.*12.*2.*M_PI);
    //reb_integrate(r, 200.);
    reb_free_simulation(r);
}

void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}
    if (reb_output_check(r, (4. / 365.25) * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");

        // rotate whole simulation to rotating frame
        struct reb_vec3d v1 = {.x = r->particles[1].x, .y = r->particles[1].y, .z = r->particles[1].z};
        struct reb_vec3d v2 = {.x = (1 / (1 + (0.01 / (1 - 0.01)))) * 5.2, .y = 0, .z = 0};
        struct reb_rotation r1 = reb_rotation_init_from_to(v1, v2);
        //reb_simulation_irotate(r, r1);
        double e = jacobi_dh(r);//reb_tools_energy(r);
        fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e\n",r->t, e, r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z);
        fclose(f);

        reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);

        struct reb_rotation inverse = reb_rotation_inverse(r1);
        //reb_simulation_irotate(r, inverse);
    }
}
