/**
 * Kozai cycles
 *
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star.
 * The integrator automatically adjusts the timestep so that
 * even very high eccentricity encounters are resolved with high
 * accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "integrator_trace.h"
#include "integrator_whfast.h"

void heartbeat(struct reb_simulation* r);

double e_init;
char title[100] = "trace_kozai";
char title_remove[100] = "rm -rf trsace_tol_kozai";
clock_t begin;
clock_t end;

double tmax = 1e4 * 2 * M_PI;
int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    //reb_simulation_start_server(r, 1234);

    // Setup constants
    //r->dt           = M_PI*1e-2;     // initial timestep
    //r->integrator   = REB_INTEGRATOR_WHFAST;
    //r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    r->integrator = REB_INTEGRATOR_TRACE;
    //r->ri_bs.eps_rel = 1e-11;            // Relative tolerance
    //r->ri_bs.eps_abs = 1e-11;            // Absolute tolerance
    //r->ri_ias15.adaptive_mode = 2;
    r->heartbeat    = heartbeat;
    r->N_active=2;

    // Initial conditions

    struct reb_particle star = {0};
    star.m  = 1;
    reb_simulation_add(r, star);

    // The planet (a zero mass test particle)
    struct reb_particle planet = {0};
    double e_testparticle = 0;
    planet.m  = 0.0;
    planet.x  = 1.-e_testparticle;
    planet.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle));

    struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);
    r->dt = o.P/100.1234;

    // The perturber
    struct reb_particle perturber = {0};
    perturber.x  = 10;
    double inc_perturber = 89.9;
    perturber.m  = 1;
    perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt((star.m+perturber.m)/perturber.x);
    perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt((star.m+perturber.m)/perturber.x);
    reb_simulation_add(r, perturber);
    reb_simulation_add(r, planet);

/*
    struct reb_particle star = {0};
    star.m = 0.809;
    star.r = 0.683*0.00465;
    reb_simulation_add(r, star);

   // HAT-P-11b
   // Yee et al 2018
   struct reb_particle planet = {0};
   double planet_m  = 0.0736 * 9.55e-4; // jupiter masses
   double planet_r = 4.36 * 4.2588e-5; // Earth radii
   double planet_a = 0.5;
   double planet_e = 0.1;
   double planet_Omega = (117.1 - 180.) * (M_PI / 180.); //reb_random_uniform(sim, 0., 2 * M_PI);
   double planet_inc = 39. * M_PI/180.;
   //double planet_inc = 106. * M_PI/180.;

   // HAT-P-11c - treated as a point particle

   //struct reb_particle perturber = {0};
   double perturber_a =  4.192;
   double perturber_e = 0.56;
   double perturber_m = 3.06 * 9.55e-4;
   double perturber_inc = 33.5 * M_PI / 180.;
   double perturber_Omega = 117.1 * M_PI / 180.;

   reb_simulation_add_fmt(r, "m r a e inc Omega", planet_m, planet_r, planet_a, planet_e, planet_inc, planet_Omega);
   reb_simulation_add_fmt(r, "m a e inc Omega", perturber_m, perturber_a, perturber_e, perturber_inc, perturber_Omega);
*/
    reb_simulation_move_to_com(r);

    remove(title_remove);        // delete previous output file
    FILE* f = fopen(title,"w");
    fprintf(f, "t,E,a1,i1,e1,a2,i2,e2,time_spent\n");
    fclose(f);

    begin = clock();
    e_init = reb_simulation_energy(r);

    while (r->t < tmax){
      FILE* f = fopen(title,"a");
      end = clock();
      double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

      struct rebx_extras* const rebx = r->extras;
      struct reb_particle* sun = &r->particles[0];
      struct reb_particle* p1 = &r->particles[1];
      struct reb_particle* pert = &r->particles[2];

      // orbits
      struct reb_orbit o1 = reb_orbit_from_particle(r->G, *p1, *sun);
      double a1 = o1.a;
      double e1 = o1.e;
      double i1 = o1.inc;
      double Om1 = o1.Omega;
      double pom1 = o1.pomega;
      struct reb_vec3d n1 = o1.hvec;

      struct reb_orbit o2 = reb_orbit_from_particle(r->G, *pert, *sun);
      double a2 = o2.a;
      double e2 = o2.e;
      double i2 = o2.inc;
      double Om2 = o2.Omega;
      double pom2 = o2.pomega;
      struct reb_vec3d n2 = o2.hvec;
      //double mi = obl(n1,n2);

      fprintf(f, "%f,%e,%f,%f,%f,%f,%f,%f,%f\n", r->t,fabs((reb_simulation_energy(r) - e_init) / e_init),a1,i1,e1,a2,i2,e2,time_spent); // print spins and orbits
      fclose(f);
      reb_simulation_integrate(r, r->t + 1e1*2*M_PI);
      //e_init = reb_simulation_energy(r);
      //printf("%0.20e\n", fabs((reb_simulation_energy(r) - e_init) / e_init));
      }

    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
    //if(reb_simulation_output_check(r, 20.*M_PI)){       // outputs to the screen
    //    reb_simulation_output_timing(r, tmax);
    //}
    /*
    if(reb_simulation_output_check(r, 12.)){            // outputs to a file
        reb_simulation_output_orbits(r, "orbits.txt");
    }
    */


}
