#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265

struct reb_vec3d reb_calculate_quad(struct reb_simulation* sim, struct reb_particle* const pi, int i, struct reb_particle* const pj, int j){
  // This calculates the quad and tidal forces on the particle pi due to particle pj. Returns the total force
  struct reb_ode* ode = *sim->odes;

  // extract relevant Parameters for p1.  SET THESE FOR NOW
  const double k1 = 0.035;// rebx_get_param(sim->extras, p1->ap, "k");

  // for p2
  const double k2 = 0.2; //rebx_get_param(sim->extras, p2->ap, "k");

  // for p3
  const double k3 = 0.2; //rebx_get_param(sim->extras, p2->ap, "k");

  // Pack these into arrays
  const double apsidals[3] = {k1, k2, k3};

  struct reb_vec3d tot_force = {0};

  if (pi != pj){

    double rj = pj->r;
    double mi = pi->m;
    double mj = pj->m;

    // Extract spin info DiffEq
    double sx = ode->y[i * 3];
    double sy = ode->y[i * 3 + 1];
    double sz = ode->y[i * 3 + 2];

    // distance vector FROM j TO i
    double dx = pj->x - pi->x;
    double dy = pj->y - pi->y;
    double dz = pj->z - pi->z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    // Unit distance vector
    double dx_hat = dx / dist;
    double dy_hat = dy / dist;
    double dz_hat = dz / dist;

    // acceleration on i due to quadrupole of j
    double quad_prefactor = (pow(rj, 5) * (1 + (mi / mj)) * apsidals[i]) / pow(dist, 4);
    double omega_dot_rhat = sx * dx_hat + sy * dy_hat + sz * dz_hat;
    double omega_squared = sx * sx + sy * sy + sz * sz;

    double qx = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dx_hat - (2 * omega_dot_rhat * sx));
    double qy = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dy_hat - (2 * omega_dot_rhat * sy));
    double qz = quad_prefactor * ((5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3))) * dz_hat - (2 * omega_dot_rhat * sz));

    // total forces
    tot_force.x = qx;
    tot_force.y = qy;
    tot_force.z = qz;
    //printf("%0.10f\t %d\t %d\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\n", sim->t, i, j, qx, tx, qy, ty, qz, tz);
    //printf("%f\t %d\t %d\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\n", sim->t, i, j, quad_prefactor, omega_dot_rhat, omega_squared, (12. * sim->G * mj / pow(dist, 3)), dx_hat, sx, sy, sz);
  }

  return tot_force;

}

struct reb_vec3d reb_calculate_tidal(struct reb_simulation* sim, struct reb_particle* const pi, int i, struct reb_particle* const pj, int j){
  // This calculates the quad and tidal forces on the particle pi due to particle pj. Returns the total force
  struct reb_ode* ode = *sim->odes;

  // extract relevant Parameters for p1.  SET THESE FOR NOW
  const double k1 = 0.035;// rebx_get_param(sim->extras, p1->ap, "k");
  const double q1 = 100000.; //rebx_get_param(sim->extras, p1->ap, "Q");

  // for p2
  const double k2 = 0.2; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q2 = 10000.; //rebx_get_param(sim->extras, p2->ap, "Q");

  // for p3
  const double k3 = 0.2; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q3 = 10000.; //rebx_get_param(sim->extras, p2->ap, "Q");

  // Pack these into arrays
  const double apsidals[3] = {k1, k2, k3};
  const double tidals[3] = {q1, q2, q3};

  struct reb_vec3d tot_force = {0};

  if (pi != pj){

    double ri = pi->r;
    double mi = pi->m;
    double mj = pj->m;

    // Orbital elements for the mutual orbit
    // This is a bit tricky b/c sometimes the star acts as the perturber: need to check against this
    double n;
    double a;
    if (j == 0) {
      struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[0], *pi);
      n = o.n;
      a = o.a;
    }
    else {
      struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[0], *pj);
      n = o.n;
      a = o.a;
    }

    // Extract spin info DiffEq
    double sx = ode->y[i * 3];
    double sy = ode->y[i * 3 + 1];
    double sz = ode->y[i * 3 + 2];

    // distance vector FROM j TO i
    double dx = pj->x - pi->x;
    double dy = pj->y - pi->y;
    double dz = pj->z - pi->z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);
    // Unit distance vector
    double dx_hat = dx / dist;
    double dy_hat = dy / dist;
    double dz_hat = dz / dist;

    // Velocity vector: j to i
    double dvx = pj->vx - pi->vx;
    double dvy = pj->vy - pi->vy;
    double dvz = pj->vz - pi->vz;

    // relevant dot and cross products
    double rhat_dot_vel = dx_hat * dvx + dy_hat * dvy + dz_hat * dvz;
    double rhat_cross_v_x = dy_hat * dvz - dz_hat * dvy;
    double rhat_cross_v_y = dz_hat * dvx - dx_hat * dvz;
    double rhat_cross_v_z = dx_hat * dvy - dy_hat * dvx;

    // first bracketed vector
    double vec1_x = 3 * rhat_dot_vel * dx_hat;
    double vec1_y = 3 * rhat_dot_vel * dy_hat;
    double vec1_z = 3 * rhat_dot_vel * dz_hat;

    // Tidal damping on i due to j
    double prefactor = -((6 * n * apsidals[i]) / tidals[i]) * (mj / mi) * pow((ri / a), 5) * pow((a / dist), 8);

    // Build vector 2 - this is the parenthesis term in (4)
    double temp_x = rhat_cross_v_x - dist * sx;
    double temp_y = rhat_cross_v_y - dist * sy;
    double temp_z = rhat_cross_v_z - dist * sz;

    double vec2_x = temp_y * dz_hat - temp_z * dy_hat;
    double vec2_y = temp_z * dx_hat - temp_x * dz_hat;
    double vec2_z = temp_x * dy_hat - temp_y * dx_hat;

    double tx = prefactor * (vec1_x + vec2_x);
    double ty = prefactor * (vec1_y + vec2_y);
    double tz = prefactor * (vec1_z + vec2_z);

    // total forces
    tot_force.x = tx;
    tot_force.y = ty;
    tot_force.z = tz;
  }

  return tot_force;

}

// LETS PROTOTYPE THIS FOR 3 BODIES
// Hacking quad and tidal forces in for the 3-body system
void additional_forces(struct reb_simulation* sim){
  // This applies the effect of all quad and tidal forces, using the helper function
  int nb = 3;

  for (int i = 0; i < nb; i++){
    for (int j = 0; j < nb; j++){
      if (i != j){
        struct reb_particle* pi = &(sim->particles[i]);
        struct reb_particle* pj = &(sim->particles[j]);

        double mi = pi->m;
        double mj = pj->m;

        struct reb_vec3d quad = reb_calculate_quad(sim, pi, i, pj, j);
        struct reb_vec3d tidal = reb_calculate_tidal(sim, pi, i, pj, j);
        
        pi->ax -= ((mj / (mi + mj)) * (quad.x + tidal.x));
        pi->ay -= ((mj / (mi + mj)) * (quad.y + tidal.y));
        pi->az -= ((mj / (mi + mj)) * (quad.z + tidal.z));
        
        pj->ax += ((mi / (mi + mj)) * (quad.x + tidal.x));
        pj->ay += ((mi / (mi + mj)) * (quad.y + tidal.y));
        pj->az += ((mi / (mi + mj)) * (quad.z + tidal.z));
      }
    }
  }

  // Migration forces
  struct reb_particle* p2 = &sim->particles[1];
  struct reb_particle* p3 = &sim->particles[2];

  if (sim->t <= 2e6 * 2 * PI){
    p2->vx *= (1 - (1e-7 * 1e-3));
    p2->vy *= (1 - (1e-7 * 1e-3));
    p2->vz *= (1 - (1e-7 * 1e-3));

    p3->vx *= (1 - (1.15 * 1e-7 * 1e-3));
    p3->vy *= (1 - (1.15 * 1e-7 * 1e-3));
    p3->vz *= (1 - (1.15 * 1e-7 * 1e-3));

  }

}


void derivatives(struct reb_ode* ode, double* const yDot, const double* const y, const double t){
    // const double omega = sqrt(k/m);
    // struct reb_orbit o = reb_tools_particle_to_orbit(ode->r->G, ode->r->particles[1], ode->r->particles[0]);
    // double forcing = sin(o.f);

    int nb = 3;
    struct reb_simulation* sim = ode->r;

    yDot[0] = 0;
    yDot[1] = 0;
    yDot[2] = 0;
    yDot[3] = 0;
    yDot[4] = 0;
    yDot[5] = 0;
    yDot[6] = 0;
    yDot[7] = 0;
    yDot[8] = 0;

    for (int i = 0; i < nb; i++){
      for (int j = 0; j < nb; j++){
          if (i != j){
            struct reb_particle* pi = &(sim->particles[i]);
            struct reb_particle* pj = &(sim->particles[j]);

            double dx = pj->x - pi->x;
            double dy = pj->y - pi->y;
            double dz = pj->z - pi->z;
            
            double mi = pi->m;
            double mj = pj->m;

            double mu_ij = -(mi * mj) / ((mi + mj));
            double moi_i = 0.4 * pi->m * pi->r * pi->r;

            struct reb_vec3d quad = reb_calculate_quad(sim, pj, j, pi, i); // Swapped for spin EOM
            struct reb_vec3d tidal = reb_calculate_tidal(sim, pi, i, pj, j);

            double tot_x = quad.x + tidal.x;
            double tot_y = quad.y + tidal.y;
            double tot_z = quad.z + tidal.z;

            yDot[i * 3] += ((mj / (mi + mj)) * ((dy * tot_z - dz * tot_y) * (mu_ij / moi_i)));
            yDot[i * 3 + 1] += ((mj / (mi + mj)) * ((dz * tot_x - dx * tot_z) * (mu_ij / moi_i)));
            yDot[i * 3 + 2] += ((mj / (mi + mj)) * ((dx * tot_y - dy * tot_x) * (mu_ij / moi_i)));
            
            yDot[j * 3] -= ((mi / (mi + mj)) * ((dy * tot_z - dz * tot_y) * (mu_ij / moi_i)));
            yDot[j * 3 + 1] -= ((mi / (mi + mj)) * ((dz * tot_x - dx * tot_z) * (mu_ij / moi_i)));
            yDot[j * 3 + 2] -= ((mi / (mi + mj)) * ((dx * tot_y - dy * tot_x) * (mu_ij / moi_i)));
          }
      }
    }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m r", 1., 0.00465);                // Central object
    reb_add_fmt(r, "m a e r inc", 5. * 3.0e-6, 0.175, 0.01, 2.5 * 4.26e-5, 1. * (PI / 180.)); // Planet 1
    reb_add_fmt(r, "m a e r inc", 5.5 * 3.0e-6, 0.235, 0.01, 2.5 * 4.26e-5, 1. * (PI / 180.));
    reb_move_to_com(r);
    r->N_active = 3;

    r->additional_forces = additional_forces;
    r->integrator = REB_INTEGRATOR_IAS15;  // Bulirsch-Stoer integrator
    //r->ri_bs.eps_rel = 1e-8;            // Relative tolerance
    //r->ri_bs.eps_abs = 1e-8;            // Absolute tolerance
    //r->dt = 1e-2;
    r->force_is_velocity_dependent = 1;

    // Here we define the spins of the three bodies
    double solar_spin_period = 20 / (2 * PI * 365);
    double solar_spin = (2 * PI) / solar_spin_period;
    double obliquity_solar = 0. * (PI / 180.);
    double res_angle_solar = 60. * (PI / 180.);
    double spin_x_solar = solar_spin * sin(obliquity_solar) * cos(res_angle_solar);
    double spin_y_solar = solar_spin * sin(obliquity_solar) * sin(res_angle_solar);
    double spin_z_solar = solar_spin * cos(obliquity_solar);

    double spin_period_1 = 5 / (365 * 2 * PI); // 5 days in reb years
    double spin_1 = (2 * PI) / spin_period_1; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_1 = 1. * (PI / 180.);
    double res_angle_1 = 120. * (PI / 180.);
    double spin_x1 = spin_1 * sin(obliquity_1) * cos(res_angle_1);
    double spin_y1 = spin_1 * sin(obliquity_1) * sin(res_angle_1);
    double spin_z1 = spin_1 * cos(obliquity_1);

    double spin_period_2 = 3 / (365 * 2 * PI); // 5 days in reb years
    double spin_2 = (2 * PI) / spin_period_2; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_2 = 1. * (PI / 180.);
    double res_angle_2 = 180. * (PI / 180.);
    double spin_x2 = spin_2 * sin(obliquity_2) * cos(res_angle_2);
    double spin_y2 = spin_2 * sin(obliquity_2) * sin(res_angle_2);
    double spin_z2 = spin_2 * cos(obliquity_2);

    struct reb_ode* spin_eom = reb_create_ode(r,9);   // Add an ODE with 2 dimensions
    spin_eom->derivatives = derivatives;              // Right hand side of the ODE
    spin_eom->y[0] = spin_x_solar;                               // Initial conditions
    spin_eom->y[1] = spin_y_solar;
    spin_eom->y[2] = spin_z_solar;

    spin_eom->y[3] = spin_x1;                               // Initial conditions
    spin_eom->y[4] = spin_y1;
    spin_eom->y[5] = spin_z1;

    spin_eom->y[6] = spin_x2;                               // Initial conditions
    spin_eom->y[7] = spin_y2;
    spin_eom->y[8] = spin_z2;

    // printf("%0.10f\t %0.10f\t %0.10f %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f \n", spin_eom->y[0], spin_eom->y[1], spin_eom->y[2], spin_eom->y[3], spin_eom->y[4], spin_eom->y[5], spin_eom->y[6], spin_eom->y[7], spin_eom->y[8]);

    r->odes = &spin_eom;

   FILE* f = fopen("out.txt","w");
    for (int i=0; i<500; i++){

        struct reb_particle sun = r->particles[0];
        struct reb_particle p1 = r->particles[1];
        struct reb_particle p2 = r->particles[2];

        struct reb_orbit o1 = reb_tools_particle_to_orbit(r->G, p1, sun);
        double a1 = o1.a;
        struct reb_orbit o2 = reb_tools_particle_to_orbit(r->G, p2, sun);
        double a2 = o2.a;

        double mag1 = sqrt(spin_eom->y[3] * spin_eom->y[3] + spin_eom->y[4] * spin_eom->y[4] + spin_eom->y[5] * spin_eom->y[5]);
        double ob1 = acos(spin_eom->y[5] / mag1) * (180 / PI);
        double mag2 = sqrt(spin_eom->y[6] * spin_eom->y[6] + spin_eom->y[7] * spin_eom->y[7] + spin_eom->y[8] * spin_eom->y[8]);
        double ob2 = acos(spin_eom->y[8] / mag2) * (180 / PI);

        //printf("torb=%.10f \t t=%.4f\t planet = %6.3f \t spin axis = %.10f %.10f %.10f \t ob = %.3f \t ang = %.10f\n", torb, r->t, o.a, dummy2.vx, dummy2.vy, dummy2.vz, obliquity, tot_angular_momentum);
        printf("t=%f\t a1 = %.6f\t a2 = %.6f\t ob1 = %.10f\t ob2 = %.10f\n", r->t / (2 * PI), a1, a2, ob1, ob2);
        //printf("t=%f\t %0.2f\t %0.2f\t %0.2f\n", r->t, spin_eom->y[3], spin_eom->y[4], spin_eom->y[5]);

        //printf("t=%f\t a1 = %.6f\t a2 = %.6f\n", r->t / (2 * PI), a1, a2);
        //printf("%f, %f, %f, %f, %f, %f\n", sun.vx, sun.vy, sun.vz, p.vx, p.vy, p.vz);
        // fprintf(f, "%.4f %.10f %.10f %.10f %.10f\n", r->t, o.a, torb, o.e, obliquity);

        reb_integrate(r,r->t+(5000 * 2 * PI));
    }
    fclose(f);
    reb_free_ode(spin_eom);
    reb_free_simulation(r);
}
