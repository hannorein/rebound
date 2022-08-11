#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double threshold = 0.0;

void compute_transformation_angles(struct reb_simulation* sim, double* theta1, double* theta2){
    // From celmech line 330
    struct reb_vec3d gtot_vec = reb_tools_angular_momentum(sim);
    double gtot = sqrt(gtot_vec.x * gtot_vec.x + gtot_vec.y * gtot_vec.y + gtot_vec.z * gtot_vec.z);
    double ghat_x = gtot_vec.x / gtot;
    double ghat_y = gtot_vec.y / gtot;
    double ghat_z = gtot_vec.z / gtot;
    double ghat_perp = sqrt(1 - ghat_z * ghat_z);
    *theta1 = M_PI / 2 - atan2(ghat_y, ghat_x);
    *theta2 = M_PI / 2 - atan2(ghat_z, ghat_perp);
}

struct reb_vec3d EulerAnglesTransform(struct reb_vec3d xyz, double Omega, double I, double omega){
    // celmech line 341
    double x = xyz.x;
    double y = xyz.y;
    double z = xyz.z;

    double s1 = sin(omega);
    double c1 = cos(omega);
    double x1 = c1 * x - s1 * y;
    double y1 = s1 * x + c1 * y;
    double z1 = z;

    double s2 = sin(I);
    double c2 = cos(I);
    double x2 = x1;
    double y2 = c2 * y1 - s2 * z1;
    double z2 = s2 * y1 + c2 * z1;

    double s3 = sin(Omega);
    double c3 = cos(Omega);
    double x3 = c3 * x2 - s3 * y2;
    double y3 = s3 * x2 + c3 * y2;
    double z3 = z2;

    struct reb_vec3d shifted = {x3, y3, z3};
    return shifted;
}

void align_simulation(struct reb_simulation* sim){
    // celmech line 360
    const int N_real = sim->N - sim->N_var;
    double theta1, theta2;
    compute_transformation_angles(sim, &theta1, &theta2);
    for (int i = 0; i <= N_real; i++){
        struct reb_particle* p = &sim->particles[i];
	struct reb_vec3d pos = {p->x, p->y, p->z};
	struct reb_vec3d vel = {p->vx, p->vy, p->vz};

	struct reb_vec3d ps = EulerAnglesTransform(pos, 0, theta2, theta1);
	struct reb_vec3d vs = EulerAnglesTransform(vel, 0, theta2, theta1);

	p->x = ps.x;
	p->y = ps.y;
	p->z = ps.z;

	p->vx = vs.x;
	p->vy = vs.y;
	p->vz = vs.z;
    }
}

double vis_viva(struct reb_simulation* sim, struct reb_particle* const pi, struct reb_particle* const source){
    // Returns semimajor axis via the vis-viva equation
  
    double dx = source->x - pi->x;
    double dy = source->y - pi->y;
    double dz = source->z - pi->z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);

    double dvx = source->vx - pi->vx;
    double dvy = source->vy - pi->vy;
    double dvz = source->vz - pi->vz;
    double vel = sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

    return ((sim->G * source->m * dist) / (2 * sim->G * source->m - dist * vel * vel));
}

struct reb_vec3d transform(double inc, double omega, struct reb_vec3d spin_inv){
    // This ts a vector from the INVARIANT frame to the PLANET frame
    double sx = spin_inv.x;
    double sy = spin_inv.y;
    double sz = spin_inv.z;

    double t[3][3];

    t[0][0] = cos(omega);
    t[0][1] = sin(omega);
    t[0][2] = 0;
    t[1][0] = -cos(inc) * sin(omega);
    t[1][1] = cos(inc) * cos(omega);
    t[1][2] = sin(inc);
    t[2][0] = sin(inc) * sin(omega);
    t[2][1] = -sin(inc) * cos(omega);
    t[2][2] = cos(inc);

    struct reb_vec3d spin_planet = {0};

    spin_planet.x = sx * t[0][0] + sy * t[0][1] + sz * t[0][2];
    spin_planet.y = sx * t[1][0] + sy * t[1][1] + sz * t[1][2];
    spin_planet.z = sx * t[2][0] + sy * t[2][1] + sz * t[2][2];

    return spin_planet;
}
struct reb_vec3d inverse_transform(double inc, double omega, struct reb_vec3d spin_planet){ 
    // This ts a vector from the PLANET frame to the INVARIANT frame
    double sx = spin_planet.x;
    double sy = spin_planet.y;
    double sz = spin_planet.z;

    double t[3][3];

    t[0][0] = cos(omega);
    t[0][1] = -cos(inc) * sin(omega);
    t[0][2] = sin(inc) * sin(omega);
    t[1][0] = sin(omega);
    t[1][1] = cos(inc) * cos(omega);
    t[1][2] = -cos(omega) * sin(inc);
    t[2][0] = 0;
    t[2][1] = sin(inc);
    t[2][2] = cos(inc);


    struct reb_vec3d spin_inv = {0};

    spin_inv.x = sx * t[0][0] + sy * t[0][1] + sz * t[0][2];
    spin_inv.y = sx * t[1][0] + sy * t[1][1] + sz * t[1][2];
    spin_inv.z = sx * t[2][0] + sy * t[2][1] + sz * t[2][2];

    return spin_inv;
}

double mean_motion(struct reb_simulation* sim, struct reb_particle* const pi, struct reb_particle* const source, double a){
    return sqrt(sim->G * (pi->m + source->m) / pow(a, 3));
}

struct reb_vec3d reb_calculate_quad(struct reb_simulation* sim, struct reb_particle* const pi, int i, struct reb_particle* const pj, int j){
  // This calculates the quad and tidal forces on the particle pi due to particle pj. Returns the total force
  // per mardling this is f^i_QD
  
  struct reb_ode* ode = *sim->odes;

  // extract relevant Parameters for p1.  SET THESE FOR NOW
  const double k1 = 0.07;// rebx_get_param(sim->extras, p1->ap, "k");

  // for p2
  const double k2 = 0.4; //rebx_get_param(sim->extras, p2->ap, "k");

  // for p3
  const double k3 = 0.4; //rebx_get_param(sim->extras, p2->ap, "k");

  // Pack these into arrays
  const double loves[3] = {k1, k2, k3};

  struct reb_vec3d tot_force = {0};

  if (pi != pj){

    double ri = pi->r;
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
    double quad_prefactor = (pow(ri, 5) * (1 + (mj / mi)) * (loves[i] / 2.)) / pow(dist, 4);
    double omega_dot_rhat = sx * dx_hat + sy * dy_hat + sz * dz_hat;
    double omega_squared = sx * sx + sy * sy + sz * sz;
    double t1 = 5 * pow(omega_dot_rhat, 2) - omega_squared - (12. * sim->G * mj / pow(dist, 3)); 
    
    double qx = quad_prefactor * (t1 * dx_hat - (2 * omega_dot_rhat * sx));
    double qy = quad_prefactor * (t1 * dy_hat - (2 * omega_dot_rhat * sy));
    double qz = quad_prefactor * (t1 * dz_hat - (2 * omega_dot_rhat * sz));
    
 //   double t2 = 2 * omega_dot_rhat * sz;
    //printf("%d %d %f\n", i, j, 12. * sim->G * mj / pow(dist, 3));
    
//    if (sim->t / (2 * PI) > threshold && j == 0 && i == 2){
//    	printf("%f,%0.15f,%0.15f,%0.15f,%0.15f\n", sim->t / (2 * PI), quad_prefactor, t1 * dz_hat, t2, dz_hat);
//	threshold += 10.0;
//    }
    
    tot_force.x = qx;
    tot_force.y = qy;
    tot_force.z = qz;
    //printf("%0.10f\t %d\t %d\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\n", sim->t, i, j, qx, tx, qy, ty, qz, tz);
    //printf("%f\t %d\t %d\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\t %0.10f\n", sim->t, i, j, quad_prefactor, omega_dot_rhat, omega_squared, (12. * sim->G * mj / pow(dist, 3)), dx_hat, sx, sy, sz);
  }

  return tot_force;

}

struct reb_vec3d reb_calculate_tidal(struct reb_simulation* sim, struct reb_particle* const pi, int i, struct reb_particle* const pj, int j){
  struct reb_vec3d temp = {0, 0, 0};
  return temp;
  // This calculates the quad and tidal forces on the particle pi due to particle pj. Returns the total force
  struct reb_ode* ode = *sim->odes;

  // extract relevant Parameters for p1.  SET THESE FOR NOW
  const double k1 = 0.07;// rebx_get_param(sim->extras, p1->ap, "k");
  const double q1 = 100000.; //rebx_get_param(sim->extras, p1->ap, "Q");

  // for p2
  const double k2 = 0.4; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q2 = 10000.; //rebx_get_param(sim->extras, p2->ap, "Q");

  // for p3
  const double k3 = 0.4; //rebx_get_param(sim->extras, p2->ap, "k");
  const double q3 = 10000.; //rebx_get_param(sim->extras, p2->ap, "Q");

  // Pack these into arrays
  const double loves[3] = {k1, k2, k3};
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
      //struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[0], *pi);
      //n = o.n;
      //a = o.a;

      a = vis_viva(sim, pi, &sim->particles[0]);
      n = mean_motion(sim, pi, &sim->particles[0], a);
    }
    else {
      //struct reb_orbit o = reb_tools_particle_to_orbit(sim->G, sim->particles[0], *pj);
      //n = o.n;
      //a = o.a;

      a = vis_viva(sim, pj, &sim->particles[0]);
      n = mean_motion(sim, pj, &sim->particles[0], a);
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
    double prefactor = -((3. * n * loves[i]) / tidals[i]) * (mj / mi) * pow((ri / a), 5) * pow((a / dist), 8);

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

        struct reb_vec3d quad = {0};
        struct reb_vec3d tidal = {0};
	
	if (j == 0){
        quad = reb_calculate_quad(sim, pi, i, pj, j);
	tidal = reb_calculate_tidal(sim, pi, i, pj, j);
        }

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
  struct reb_particle* star = &sim->particles[0];
  struct reb_particle* p2 = &sim->particles[1];
  struct reb_particle* p3 = &sim->particles[2];

  double tau_2 = 5e6 * 2 * M_PI;
  double tau_3 = tau_2 / 1.1;

  if (sim->t <= 2e6 * 2 * M_PI){
    double vx2 = p2->vx - star->vx;
    double vy2 = p2->vy - star->vy;
    double vz2 = p2->vz - star->vz;

    double vx3 = p3->vx - star->vx;
    double vy3 = p3->vy - star->vy;
    double vz3 = p3->vz - star->vz;
    
    p2->ax -= (vx2 / (2 * tau_2)); 
    p2->ay -= (vy2 / (2 * tau_2)); 
    p2->az -= (vz2 / (2 * tau_2)); 

    p3->ax -= (vx3 / (2 * tau_3)); 
    p3->ay -= (vy3 / (2 * tau_3)); 
    p3->az -= (vz3 / (2 * tau_3)); 

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
    yDot[9] = 0;
    yDot[10] = 0;

    for (int i = 0; i < nb; i++){
      for (int j = 0; j < nb; j++){
          if (i != j){
            // all this lookup can be made more efficient
	    struct reb_particle* pi = &(sim->particles[i]);
            struct reb_particle* pj = &(sim->particles[j]);

            double dx = pj->x - pi->x;
            double dy = pj->y - pi->y;
            double dz = pj->z - pi->z;
            
            double mi = pi->m;
            double mj = pj->m;
            
            double mu_ij = -(mi * mj) / ((mi + mj));
            double moi_i;

            if (i == 0){ 
                moi_i = 0.07 * pi->m * pi->r * pi->r;
            }
            else {
                moi_i = 0.25 * pi->m * pi->r * pi->r;
            }

            struct reb_vec3d quad = reb_calculate_quad(sim, pi, i, pj, j); // Swapped for spin EOM
            struct reb_vec3d tidal = {0};
//
	    if (j == 0){
	        tidal = reb_calculate_tidal(sim, pi, i, pj, j);
	    }

            double tot_x = quad.x + tidal.x;
            double tot_y = quad.y + tidal.y;
            double tot_z = quad.z + tidal.z;
            //printf("%d %d %0.15f %0.15f\n", i, j, temp * pow(2 * PI, 2), moi_i / pow(2 * PI, 1.333333));

            yDot[i * 3] += ((dy * tot_z - dz * tot_y) * (mu_ij / moi_i));
            yDot[i * 3 + 1] += ((dz * tot_x - dx * tot_z) * (mu_ij / moi_i));
	    yDot[i * 3 + 2] += ((dx * tot_y - dy * tot_x) * (mu_ij / moi_i));
          }
      }
      // Track spin of the planet
 //     if (i != 0){
//	int index = 2 - i;
//	double result = (sqrt(ode->y[i*3] * ode->y[i*3] + ode->y[i*3+1] * ode->y[i*3+1] + ode->y[i*3+2] * ode->y[i*3+2]));
//	if (result > 2. * PI){
//	    result = fmod(result, 2. * PI);
//	}

//	yDot[10 - index] = result;
//      }
    }

 //   if (sim->t / (2 * PI) > threshold){
 //   	printf("%f,%0.15f,%0.15f,%0.15f\n", sim->t / (2 * PI), yDot[6], yDot[7], yDot[8]);
//	threshold += 10.0;
//    }
    
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m r", 1., 0.00465);                // Central object
    reb_add_fmt(r, "m a e r inc Omega pomega l", 5. * 3.0e-6, 0.173087, 0.01, 2.5 * 4.26e-5, 0.5 * (M_PI / 180.), 0.0 * (M_PI / 180.), 0.0 * (M_PI / 180.), 45.0 * (M_PI / 180.)); // Planet 1
    double p2_inc = -0.431;
    double p2_Omega = 0.0;
    reb_add_fmt(r, "m a e r inc Omega pomega l", 5. * 3.0e-6, 0.2329, 0.01, 2.5 * 4.26e-5, p2_inc * (M_PI / 180.), p2_Omega * (M_PI / 180.), 0.0 * (M_PI / 180.), 90.0 * (M_PI / 180.));
    reb_move_to_com(r);
    r->N_active = 3;
    //align_simulation(r);

    //struct reb_orbit o2_0 = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);
    //p2_inc = o2_0.inc;
    //p2_Omega = o2_0.Omega;
    //printf("%f %f\n", p2_inc * (M_PI / 180), p2_Omega * (M_PI / 180));

    r->additional_forces = additional_forces;
    r->integrator = REB_INTEGRATOR_WHFAST;  // Bulirsch-Stoer integrator
  //  r->ri_bs.eps_rel = 1e-8;            // Relative tolerance
  //  r->ri_bs.eps_abs = 1e-8;            // Absolute tolerance
//    r->ri_whfast.safe_mode = 0; 
    r->dt = 1e-2;
    r->force_is_velocity_dependent = 1;
    

    // Here we define the spins of the three bodies
    struct reb_vec3d solar_sv = {0};
    double solar_spin_period = 20 * 2 * M_PI / 365;
    double solar_spin = (2 * M_PI) / solar_spin_period;
    double obliquity_solar = 0. * (M_PI / 180.);
    double res_angle_solar = 60. * (M_PI / 180.);
    double solar_c = 1;
    solar_sv.x = solar_c * solar_spin * 0.0; //sin(obliquity_solar) * cos(res_angle_solar);
    solar_sv.y = solar_c * solar_spin * 0.0; //sin(obliquity_solar) * sin(res_angle_solar);
    solar_sv.z = solar_c * solar_spin * 1.0; //cos(obliquity_solar);
//    printf("%f\n", solar_spin);
    //printf("%f %f %f\n", solar_sv.x, solar_sv.y, solar_sv.z);
    
    double planet_c = 1;
    struct reb_vec3d p1_sv_p = {0};
    double spin_period_1 = 5 * 2 * M_PI / 365; // 5 days in reb years
    double spin_1 = (2 * M_PI) / spin_period_1; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_1 = 1. * (M_PI / 180.);
    double res_angle_1 = 120. * (M_PI / 180.);
    p1_sv_p.x = planet_c * spin_1 * 0.0; //sin(obliquity_1) * cos(res_angle_1);
    p1_sv_p.y = planet_c * spin_1 * -0.0261769; //sin(obliquity_1) * sin(res_angle_1);
    p1_sv_p.z = planet_c * spin_1 * 0.99965732; // cos(obliquity_1);
    struct reb_vec3d p1_spin = {p1_sv_p.x, p1_sv_p.y, p1_sv_p.z}; // inverse_transform(1. * (PI / 180.), 40. * (PI / 180), p1_sv_p);
//    printf("%f %f %f\n", p1_spin.x * 2 * PI, p1_spin.y * 2 * PI, p1_spin.z * 2 * PI);
    

    struct reb_vec3d p2_sv_p = {0};
    double spin_period_2 = 3 * 2 * M_PI / 365; // 5 days in reb years
    double spin_2 = (2 * M_PI) / spin_period_2; // magnitude of spin vector is spin PERIOD in rad/reb years
    double obliquity_2 = 15. * (M_PI / 180.);
    double res_angle_2 = 90. * (M_PI / 180.);
    p2_sv_p.x = spin_2 * 0.0;//* sin(obliquity_2) * cos(res_angle_2);
    p2_sv_p.y = spin_2 * 0.0249736362618; //* sin(obliquity_2) * sin(res_angle_2);
    p2_sv_p.z = spin_2 * 0.999688110108; //* cos(obliquity_2);
    struct reb_vec3d p2_spin = {p1_sv_p.x, p1_sv_p.y, p1_sv_p.z}; //inverse_transform(p2_inc, p2_Omega, p2_sv_p);
//    printf("%f %f %f\n", p2_spin.x * 2 * PI, p2_spin.y * 2 * PI, p2_spin.z * 2 * PI);
    
    struct reb_ode* spin_eom = reb_create_ode(r,9);   // Add an ODE with 9 dimensions. 3 for each body's spin axis
    spin_eom->derivatives = derivatives;              // Right hand side of the ODE
    spin_eom->y[0] = solar_sv.x;                               // Initial conditions
    spin_eom->y[1] = solar_sv.y;
    spin_eom->y[2] = solar_sv.z;

    spin_eom->y[3] = p1_spin.x;                               // Initial conditions
    spin_eom->y[4] = p1_spin.y;
    spin_eom->y[5] = p1_spin.z;

    spin_eom->y[6] = p2_spin.x;                               // Initial conditions
    spin_eom->y[7] = p2_spin.y;
    spin_eom->y[8] = p2_spin.z;

    //spin_eom->y[9] = 0.0; // Initial spin configuration psi = 0
    //spin_eom->y[10] = 0.0;

    r->odes = &spin_eom;
    //clock_t start, end;
    //double cpu_time_used;

   FILE* f = fopen("8_11_sm_test.txt","w");
   fprintf(f, "t,a1,i1,e1,sx1,sy1,sz1,S1,omega1,Omega1,f1,x1,y1,z1,a2,i2,e2,sx2,sy2,sz2,S2,omega2,Omega2,f2,x2,y2,z2\n");
   //fprintf(f, "t,a1,a2\n");

//   FILE* star_f1 = fopen("4_21_forces_star_on_p1_no_dsz.txt", "w");
//   fprintf(star_f1, "t,qx,qy,qz,tx,ty,tz,x,y,z\n");
   
//   FILE* star_f2 = fopen("4_21_forces_star_on_p2_no_dsz.txt", "w");
//   fprintf(star_f2, "t,qx,qy,qz,tx,ty,tz,x,y,z\n");
   
//   FILE* p1_f = fopen("4_21_forces_p2_on_p1_no_dsz.txt", "w");
//   fprintf(p1_f, "t,qx,qy,qz,x,y,z\n");
   
//   FILE* p2_f = fopen("4_21_forces_p1_on_p2_no_dsz.txt", "w");
//   fprintf(p2_f, "t,qx,qy,qz,x,y,z\n");

    //start = clock();
    for (int i=0; i<100000; i++){
	
        struct reb_particle sun = r->particles[0];
        struct reb_particle p1 = r->particles[1];
        struct reb_particle p2 = r->particles[2];

        struct reb_orbit o1 = reb_tools_particle_to_orbit(r->G, p1, sun);
        double a1 = o1.a;//vis_viva(r, &p1, &sun);
	double O1 = o1.Omega;
	double i1 = o1.inc;
	double pom1 = o1.pomega;
	double f1 = o1.f;
	double e1 = o1.e;

        struct reb_orbit o2 = reb_tools_particle_to_orbit(r->G, p2, sun);
        double a2 = o2.a;//vis_viva(r, &p2, &sun);
	double O2 = o2.Omega;
	double i2 = o2.inc;
	double pom2 = o2.pomega;
	double f2 = o2.f;
	double e2 = o2.e;

        struct reb_vec3d s1_inv = {spin_eom->y[3], spin_eom->y[4], spin_eom->y[5]};
        struct reb_vec3d s2_inv = {spin_eom->y[6], spin_eom->y[7], spin_eom->y[8]};

        struct reb_vec3d s1 = transform(i1, O1, s1_inv);
        struct reb_vec3d s2 = transform(i2, O2, s2_inv);

        // Interpret in the planet frame
        double mag1 = sqrt(s1.x * s1.x + s1.y * s1.y + s1.z * s1.z);
        double ob1 = acos(s1.z / mag1) * (180 / M_PI);
        double mag2 = sqrt(s2.x * s2.x + s2.y * s2.y + s2.z * s2.z);
        double ob2 = acos(s2.z / mag2) * (180 / M_PI);

	struct reb_vec3d l = reb_tools_angular_momentum(r);
	double l_mag = sqrt(l.x * l.x + l.y * l.y + l.z * l.z);
        
//	double dx_1 = sun.x - p1.x;
//        double dy_1 = sun.y - p1.y;
//        double dz_1 = sun.z - p1.z;
	
//	double dx_2 = sun.x - p2.x;
//        double dy_2 = sun.y - p2.y;
//        double dz_2 = sun.z - p2.z;

//	double dx_pp = p1.x - p2.x;
//        double dy_pp = p1.y - p2.y;
//        double dz_pp = p1.z - p2.z;

        
//        double a1_v = vis_viva(r, &p1, &sun);
//        double a2_v = vis_viva(r, &p2, &sun);
         
        if (i % 1000 == 0){
            printf("t=%f\t a1=%.6f\t a2=%.6f\t o1=%0.5f\t o2=%0.5f, L = %f %f %f\n", r->t / (2 * M_PI), a1, a2, ob1, ob2, l.x / l_mag, l.y / l_mag, l.z / l_mag);
	    //printf("%f\n", r->t / (2 * PI));
        }
        //printf("t = %f\t a1 = %.10f\t a1v = %0.10f\t dif_1 = %0.5f\t a2 = %0.10f\t a2v = %0.10f\t dif_2 = %0.5f\n", r->t / (2 * PI), a1, a1_v, (a1_v - a1) / a1 * 100.0, a2, a2_v, (a2_v - a2) / a2 * 100.0);
        fprintf(f, "%f,%f,%f,%f,%0.10f,%0.10f,%0.10f,%0.10f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%0.10f,%0.10f,%0.10f,%0.10f,%f,%f,%f,%f,%f,%f\n", r->t / (2 * M_PI), a1, i1, e1, s1.x, s1.y, s1.z, mag1, pom1, O1, f1, p1.x,p1.y,p1.z,a2, i2, e2, s2.x, s2.y, s2.z, mag2, pom2, O2, f2, p2.x,p2.y,p2.z);
	//fprintf(f, "%f,%f,%f\n", r->t / (2 * PI), a1, a2);
        //printf("t=%f\t %0.2f\t %0.2f\t %0.2f\n", r->t, spin_eom->y[3], spin_eom->y[4], spin_eom->y[5]);

  //      struct reb_vec3d quad_s_2 = reb_calculate_quad(r, &p2, 2, &sun, 0);
//	struct reb_vec3d tidal_s_2 = reb_calculate_tidal(r, &p2, 2, &sun, 0);
        
//	struct reb_vec3d quad_s_1 = reb_calculate_quad(r, &p1, 1, &sun, 0);
//	struct reb_vec3d tidal_s_1 = reb_calculate_tidal(r, &p1, 1, &sun, 0);
        
//	struct reb_vec3d quad_p1_on_p2 = reb_calculate_quad(r, &p2, 2, &p1, 1);
//	struct reb_vec3d quad_p2_on_p1 = reb_calculate_quad(r, &p1, 1, &p2, 2);

	//fprintf(star_f1, "%f,%0.15f,%0.15f,%0.15f,%0.15f,%0.15f,%0.15f,%f,%f,%f\n", r->t / (2 * PI), quad_s_1.x, quad_s_1.y, quad_s_1.z, tidal_s_1.x, tidal_s_1.y, tidal_s_1.z, dx_1, dy_1, dz_1);
	//fprintf(star_f2, "%f,%0.15f,%0.15f,%0.15f,%0.15f,%0.15f,%0.15f,%f,%f,%f\n", r->t / (2 * PI), quad_s_2.x, quad_s_2.y, quad_s_2.z, tidal_s_2.x, tidal_s_2.y, tidal_s_2.z, dx_2, dy_2, dz_2);
	//fprintf(p1_f, "%f,%0.15f,%0.15f,%0.15f,%f,%f,%f\n", r->t / (2 * PI), quad_p2_on_p1.x, quad_p2_on_p1.y, quad_p2_on_p1.z, -dx_pp, -dy_pp, -dz_pp);
	//fprintf(p2_f, "%f,%0.15f,%0.15f,%0.15f,%f,%f,%f\n", r->t / (2 * PI), quad_p1_on_p2.x, quad_p1_on_p2.y, quad_p1_on_p2.z, dx_pp, dy_pp, dz_pp);
        //fprintf(f, "%.4f %.10f %.10f %.10f %.10f\n", r->t, o.a, torb, o.e, obliquity);
        reb_integrate(r,r->t+(40 * 2 * M_PI));
    }

    fclose(f);
    reb_free_ode(spin_eom);
    reb_free_simulation(r);
    printf("SUCCESS\n");
}
