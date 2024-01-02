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
 #include "integrator_trace.h"

void heartbeat(struct reb_simulation* r);
double tmax = 5e7 * 2 * M_PI;
clock_t begin;
clock_t end;
double e_init;
char title[100] = "trace_nt";
char title_remove[100] = "rm -rf trace_nt";

double obl(struct reb_vec3d v1, struct reb_vec3d v2){
  return acos(reb_vec3d_dot(v1,v2) / (sqrt(reb_vec3d_length_squared(v1)) * sqrt(reb_vec3d_length_squared(v2))));
}

int main(int argc, char* argv[]){
    // Setup constants
    struct reb_simulation* sim = reb_simulation_create();
    //sim->dt             = 2*M_PI*1.;     // initial timestep
    sim->integrator        = REB_INTEGRATOR_TRACE;
    //sim->ri_ias15.adaptive_mode=2;
    sim->heartbeat        = heartbeat;

    // Initial conditions Naoz et al 2016 Figure 16
    struct reb_particle star = {0};
    star.m  = 0.32;
    star.r = 0.1 * 0.00465;
    reb_simulation_add(sim, star);

    double planet_m  = 5.15e-5; // in Jupiter masses
    double planet_r = 1.65e-4;//4.676e-4;
    double planet_a = 2.;
    double planet_e = 0.01;
    double planet_omega = 0.;
    reb_simulation_add_fmt(sim, "m r a e", planet_m, planet_r, planet_a, planet_e);

    struct reb_orbit o = reb_orbit_from_particle(sim->G, sim->particles[1], sim->particles[0]);
    sim->dt = o.P/50.12345;
    sim->ri_trace.S_peri = reb_integrator_trace_peri_switch_default;
    sim->ri_trace.peri_distance = 1.;
    //sim->ri_trace.pfdot=1000.;
    //sim->ri_bs.eps_rel = 1e-11;            // Relative tolerance
    //sim->ri_bs.eps_abs = 1e-11;            // Absolute tolerance
    sim->exact_finish_time=0;

    // The perturber
    double perturber_mass = 10. * 9.55e-4;
    double perturber_a  = 50.;
    double perturber_e = 0.52;
    double perturber_inc = 65. * M_PI/180.;
    double perturber_omega = 0.;

    reb_simulation_add_fmt(sim, "m a e inc omega", perturber_mass, perturber_a, perturber_e, perturber_inc, perturber_omega);

    struct reb_vec3d newz = reb_simulation_angular_momentum(sim);
    struct reb_vec3d newx = reb_vec3d_cross((struct reb_vec3d){.z =1}, newz);
    struct reb_rotation rot = reb_rotation_init_to_new_axes(newz, newx);
    reb_simulation_irotate(sim, rot);
    reb_simulation_move_to_com(sim);

    // AND THIS FOR TIDES_SPIN

    system(title_remove);
    FILE* f = fopen(title,"w");
    //fprintf(f, "t,E,a1,i1,e1,p_ob,a2,i2,e2,pert_ob,mi,theta_p,time_spent\n");
    fprintf(f, "t,E,a1,i1,e1,a2,i2,e2,mi,time_spent,hj,hk\n");
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

      reb_integrator_trace_inertial_to_dh(sim);
      double px0=0., py0=0., pz0=0.;
      double mtot = sim->particles[0].m;

      for (unsigned int i=0;i<sim->N;i++){
          px0 += sim->particles[i].vx*sim->particles[i].m;
          py0 += sim->particles[i].vy*sim->particles[i].m;
          pz0 += sim->particles[i].vz*sim->particles[i].m;
          mtot += sim->particles[i].m;
      }

      int j = 1;
      const double pjx = (sim->particles[j].vx * sim->particles[j].m) - (sim->particles[j].m / mtot) * px0;
      const double pjy = (sim->particles[j].vy * sim->particles[j].m) - (sim->particles[j].m / mtot) * py0;
      const double pjz = (sim->particles[j].vz * sim->particles[j].m) - (sim->particles[j].m / mtot) * pz0;
      const double p02 = px0*px0 + py0*py0+pz0*pz0;
      const double pj2 = pjx*pjx + pjy * pjy + pjz * pjz;
      const double hj = p02 / (2. * sim->particles[0].m);

      // Kepler
      const double dx = sim->particles[j].x - sim->particles[0].x;
      const double dy = sim->particles[j].y - sim->particles[0].y;
      const double dz = sim->particles[j].z - sim->particles[0].z;
      const double d2 = dx*dx + dy*dy + dz*dz;
      const double d = sqrt(d2);
      const double hk = pj2 / (2 * sim->particles[j].m) - ((sim->G * sim->particles[0].m * sim->particles[j].m) / d);

      fprintf(f, "%f,%e,%f,%f,%f,%f,%f,%f,%f,%f,%e,%e\n", sim->t,fabs((reb_simulation_energy(sim) - e_init) / e_init),a1,i1,e1,a2,i2,e2,mi,time_spent,hj,hk); // print spins and orbits
      fclose(f);
      reb_integrator_trace_dh_to_inertial(sim);
      reb_simulation_integrate(sim, sim->t + 1e4*2*M_PI);
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