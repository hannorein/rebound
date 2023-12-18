/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 *
 * This example includes self-consistent spin, tidal & dynamical effects
 * as well as general relativity
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include <unistd.h>
 #include <math.h>
 #include <time.h>
 #include "rebound.h"

void heartbeat(struct reb_simulation* r);
double tmax = 1e7 * 2 * M_PI;
clock_t begin;
clock_t end;
double e_init;
char title[100] = "ias15_kozai";
char title_remove[100] = "rm -rf ias15_kozai";

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_simulation_create();
    // Setup constants
    //sim->dt             = 2*M_PI*1.;     // initial timestep
    sim->integrator        = REB_INTEGRATOR_IAS15;
    //sim->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    sim->ri_ias15.adaptive_mode=2;
    sim->heartbeat        = heartbeat;

    // Initial conditions

    struct reb_particle star = {0};
    star.m  = 1.1;
    star.r = 0.00465;
    reb_simulation_add(sim, star);

    double planet_m  = 7.8 * 9.55e-4; // in Jupiter masses
    double planet_r = 4.676e-4;
    double planet_a = 5.0;
    double planet_e = 0.1;
    reb_simulation_add_fmt(sim, "m r a e", planet_m, planet_r, planet_a, planet_e);

    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = o.P/100.12345;
    sim->exact_finish_time=0;

    // The perturber
    double perturber_inc = 89.6 * (M_PI / 180.);
    double perturber_mass = 0.5;
    double perturber_a  = 500.;
    double perturber_e = 0.5;
    reb_simulation_add_fmt(sim, "m a e inc", perturber_mass, perturber_a, perturber_e, perturber_inc);
    reb_simulation_move_to_com(sim);

    struct reb_vec3d newz = reb_simulation_angular_momentum(sim);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    reb_simulation_irotate(sim, rot);

    // AND THIS FOR TIDES_SPIN

    system(title_remove);
    FILE* f = fopen(title,"w");
    //fprintf(f, "t,E,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,theta_p,time_spent\n");
    fprintf(f, "t,E,a1,i1,e1,a2,i2,e2,mi,time_spent\n");
    fclose(f);

    begin = clock();
    e_init = reb_simulation_energy(sim);

    while (sim->t < tmax){
      FILE* f = fopen(title,"a");
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

      struct reb_particle* sun = &sim->particles[0];
      struct reb_particle* p1 = &sim->particles[1];
      struct reb_particle* pert = &sim->particles[2];

      // orbits
      struct reb_orbit o1 = reb_orbit_from_particle(sim->G, *p1, *sun);
      double a1 = o1.a;
      double e1 = o1.e;
      double i1 = o1.inc;
      double Om1 = o1.Omega;
      double pom1 = o1.pomega;
      struct reb_vec3d n1 = o1.hvec;

      struct reb_orbit o2 = reb_orbit_from_particle(sim->G, *pert, *sun);
      double a2 = o2.a;
      double e2 = o2.e;
      double i2 = o2.inc;
      double Om2 = o2.Omega;
      double pom2 = o2.pomega;
      struct reb_vec3d n2 = o2.hvec;
      double mi = obl(n1,n2);

      fprintf(f, "%f,%e,%f,%f,%f,%f,%f,%f,%f,%f\n", sim->t,fabs((reb_simulation_energy(sim) - e_init) / e_init),a1,i1,e1,a2,i2,e2,mi,time_spent); // print spins and orbits
      fclose(f);
      reb_simulation_integrate(sim, sim->t + 1e5*2*M_PI);
      }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    reb_simulation_free(sim);
}

void heartbeat(struct reb_simulation* sim){
   // Output spin and orbital information to file

   if (reb_simulation_output_check(sim, 1e3*2.*M_PI)){
        reb_simulation_output_timing(sim, tmax);
   }
}
