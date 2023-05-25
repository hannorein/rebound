/**
 * @file    integrator_mercurius.c
 * @brief   MERCURIUS, a modified version of John Chambers' MERCURY algorithm
 *          using the IAS15 integrator and WHFast. It works with planet-planry
 *          collisions, test particles, and additional forces.
 * @author  Hanno Rein, Dan Tamayo
 *
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein, Dan Tamayo
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_mercurius.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "integrator_bs.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

double reb_integrator_mercurius_F_cond(double d, double dcrit){
  return d - dcrit;
}
double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function used by the Mercury integrator.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return 10.*(y*y*y) - 15.*(y*y*y*y) + 6.*(y*y*y*y*y);
    }
}

double reb_integrator_mercurius_L_C4(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function C4 proposed by Hernandez (2019)
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return (70.*y*y*y*y -315.*y*y*y +540.*y*y -420.*y +126.)*y*y*y*y*y;
    }
}

double reb_integrator_mercurius_L_C5(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function C5 proposed by Hernandez (2019)
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return (-252.*y*y*y*y*y +1386.*y*y*y*y -3080.*y*y*y +3465.*y*y -1980.*y +462.)*y*y*y*y*y*y;
    }
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}

double reb_integrator_mercurius_L_infinity(const struct reb_simulation* const r, double d, double dcrit){
    // Infinitely differentiable function.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
}

// TLu F condition from Hernandez & Dehnen 2023
void reb_integrator_mercurius_inertial_to_dh(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_vec3d com_pos = {0};
    struct reb_vec3d com_vel = {0};
    double mtot = 0.;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;
    const int N = r->N;
    for (int i=0;i<N_active;i++){
        double m = particles[i].m;
        com_pos.x += m * particles[i].x;
        com_pos.y += m * particles[i].y;
        com_pos.z += m * particles[i].z;
        com_vel.x += m * particles[i].vx;
        com_vel.y += m * particles[i].vy;
        com_vel.z += m * particles[i].vz;
        mtot += m;
    }
    com_pos.x /= mtot; com_pos.y /= mtot; com_pos.z /= mtot;
    com_vel.x /= mtot; com_vel.y /= mtot; com_vel.z /= mtot;
    // Particle 0 is also changed to allow for easy collision detection
    for (int i=N-1;i>=0;i--){
        particles[i].x -= particles[0].x;
        particles[i].y -= particles[0].y;
        particles[i].z -= particles[0].z;
        particles[i].vx -= com_vel.x;
        particles[i].vy -= com_vel.y;
        particles[i].vz -= com_vel.z;
    }
    r->ri_mercurius.com_pos = com_pos;
    r->ri_mercurius.com_vel = com_vel;
}

void reb_integrator_mercurius_dh_to_inertial(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle temp = {0};
    const int N = r->N;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;
    for (int i=1;i<N_active;i++){
        double m = particles[i].m;
        temp.x += m * particles[i].x;
        temp.y += m * particles[i].y;
        temp.z += m * particles[i].z;
        temp.vx += m * particles[i].vx;
        temp.vy += m * particles[i].vy;
        temp.vz += m * particles[i].vz;
        temp.m += m;
    }
    temp.m += r->particles[0].m;
    temp.x /= temp.m;
    temp.y /= temp.m;
    temp.z /= temp.m;
    temp.vx /= particles[0].m;
    temp.vy /= particles[0].m;
    temp.vz /= particles[0].m;
    // Use com to calculate central object's position.
    // This ignores previous values stored in particles[0].
    // Should not matter unless collisions occured.
    particles[0].x = r->ri_mercurius.com_pos.x - temp.x;
    particles[0].y = r->ri_mercurius.com_pos.y - temp.y;
    particles[0].z = r->ri_mercurius.com_pos.z - temp.z;

    for (int i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += r->ri_mercurius.com_vel.x;
        particles[i].vy += r->ri_mercurius.com_vel.y;
        particles[i].vz += r->ri_mercurius.com_vel.z;
    }
    particles[0].vx = r->ri_mercurius.com_vel.x - temp.vx;
    particles[0].vy = r->ri_mercurius.com_vel.y - temp.vy;
    particles[0].vz = r->ri_mercurius.com_vel.z - temp.vz;
}


static void reb_mercurius_encounter_predict(struct reb_simulation* const r){
  // For now this flags ALL particles TLu
    struct reb_simulation_integrator_mercurius* rim = &(r->ri_mercurius);
    struct reb_particle* const particles = r->particles;
    struct reb_particle* const particles_backup = rim->particles_backup;
    const double* const dcrit = rim->dcrit;
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const double dt = r->dt;
    rim->encounterN = 1;
    rim->encounter_map[0] = 1;
    if (r->testparticle_type==1){
        rim->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        rim->tponly_encounter = 1;
    }
    for (int i=1; i<N; i++){
        rim->encounter_map[i] = 0;
    }
    for (int i=0; i<N_active; i++){
        for (int j=i+1; j<N; j++){
            if (i != 0){ // if encounter_predict is called we need to integrate ALL particles
                if (rim->encounter_map[i]==0){
                    rim->encounter_map[i] = i;
                    rim->encounterN++;
                }
                if (rim->encounter_map[j]==0){
                    rim->encounter_map[j] = j;
                    rim->encounterN++;

                }
                if (j<N_active){ // Two massive particles have a close encounter
                    rim->tponly_encounter = 0;
                }
            }
        }
    }
}

void reb_integrator_mercurius_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_mercurius_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;

    struct reb_simulation_integrator_mercurius* rim = &(r->ri_mercurius);

    const int N_active = r->N_active==-1?r->N:r->N_active;
    const int N = r->testparticle_type==0 ? N_active: r->N;
    double px=0., py=0., pz=0.;
    for (int i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m;
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px /= r->particles[0].m;
    py /= r->particles[0].m;
    pz /= r->particles[0].m;
    const int N_all = r->N;
    int current_L;
    for (int i=1;i<N_all;i++){
        current_L = rim->current_Ls[i-1];
        particles[i].x += dt*px * (1 - current_L * (particles[i].m!=0)); // TLu crude test particle fix
        particles[i].y += dt*py * (1 - current_L * (particles[i].m!=0));
        particles[i].z += dt*pz * (1 - current_L * (particles[i].m!=0));
    }
}

void reb_integrator_mercurius_com_step(struct reb_simulation* const r, double dt){
    //printf("COM\n");
    r->ri_mercurius.com_pos.x += dt*r->ri_mercurius.com_vel.x;
    r->ri_mercurius.com_pos.y += dt*r->ri_mercurius.com_vel.y;
    r->ri_mercurius.com_pos.z += dt*r->ri_mercurius.com_vel.z;
}

// Old Kepler
void reb_integrator_mercurius_kepler_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        reb_whfast_kepler_solver(r,r->particles,r->G*r->particles[0].m,i,dt); // in dh
    }
}


static void reb_mercurius_bs_encounter_step(struct reb_simulation* const r, const double _dt){
  //printf("BS Encounter\n");
  struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
  rim->mode = 1;
  // run
  const double old_dt = r->dt;
  const double old_t = r->t;
  double t_needed = r->t + _dt;
  reb_integrator_bs_reset(r);

  r->dt = 0.0001*_dt; // start with a small timestep.
  while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 ){
      struct reb_particle star = r->particles[0]; // backup velocity
      r->particles[0].vx = 0; // star does not move in dh
      r->particles[0].vy = 0;
      r->particles[0].vz = 0;
      reb_update_acceleration(r);
      reb_integrator_bs_part2(r);

      r->particles[0].vx = star.vx; // restore every timestep for collisions
      r->particles[0].vy = star.vy;
      r->particles[0].vz = star.vz;

      if (r->t+r->dt >  t_needed){
          r->dt = t_needed-r->t;
      }

      // Search and resolve collisions
      reb_collision_search(r);

      // Do any additional post_timestep_modifications.
      // Note: post_timestep_modifications is called here but also
      // at the end of the full timestep. The function thus needs
      // to be implemented with care as not to do the same
      // modification multiple times. To do that, check the value of
      // r->ri_mercurius.mode
      if (r->post_timestep_modifications){
          r->post_timestep_modifications(r);
      }

      star.vx = r->particles[0].vx; // keep track of changed star velocity for later collisions
      star.vy = r->particles[0].vy;
      star.vz = r->particles[0].vz;
      if (r->particles[0].x !=0 || r->particles[0].y !=0 || r->particles[0].z !=0){
          // Collision with star occured
          // Shift all particles back to heliocentric coordinates
          // Ignore stars velocity:
          //   - will not be used after this
          //   - com velocity is unchained. this velocity will be used
          //     to reconstruct star's velocity later.
          for (int i=r->N-1; i>=0; i--){
              r->particles[i].x -= r->particles[0].x;
              r->particles[i].y -= r->particles[0].y;
              r->particles[i].z -= r->particles[0].z;
          }
      }
    }

    if(rim->tponly_encounter){
        for (int i=1;i<rim->encounterNactive;i++){
            unsigned int mi = i;//rim->encounter_map[i];
            r->particles[mi] = rim->particles_backup[mi];
        }
    }

    r->t = old_t;
    r->dt = old_dt;
    rim->mode = 0;
}

double reb_integrator_mercurius_calculate_dcrit_for_particle(struct reb_simulation* r, unsigned int i){
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    const double m0 = r->particles[0].m;
    const double dx  = r->particles[i].x;  // in dh
    const double dy  = r->particles[i].y;
    const double dz  = r->particles[i].z;
    const double dvx = r->particles[i].vx - r->particles[0].vx;
    const double dvy = r->particles[i].vy - r->particles[0].vy;
    const double dvz = r->particles[i].vz - r->particles[0].vz;
    const double _r = sqrt(dx*dx + dy*dy + dz*dz);
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;

    const double GM = r->G*(m0+r->particles[i].m);
    const double a = GM*_r / (2.*GM - _r*v2);
    const double vc = sqrt(GM/fabs(a));
    double dcrit = 0;
    // Criteria 1: average velocity
    dcrit = MAX(dcrit, vc*0.4*r->dt);
    // Criteria 2: current velocity
    dcrit = MAX(dcrit, sqrt(v2)*0.4*r->dt);
    // Criteria 3: Hill radius
    dcrit = MAX(dcrit, rim->hillfac*a*cbrt(r->particles[i].m/(3.*r->particles[0].m)));
    // Criteria 4: physical radius
    dcrit = MAX(dcrit, 2.*r->particles[i].r);
    return dcrit;
}


void reb_integrator_mercurius_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }

    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    const int N = r->N;
    //rim->encounterN = N; // TLu temp for full IAS15 integration

    if (rim->dcrit_allocatedN<N){
        // Need to safe these arrays in SimulationArchive
        rim->dcrit              = realloc(rim->dcrit, sizeof(double)*N);
        rim->dcrit_allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_dcrit_this_timestep        = 1;
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
        rim->recalculate_coordinates_this_timestep = 1;
    }
    if (rim->allocatedN<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        rim->particles_backup       = realloc(rim->particles_backup,sizeof(struct reb_particle)*N);

        rim->current_Ks = realloc(rim->current_Ks, sizeof(int)*((N-1)*(N-2))/2);
        rim->current_Ls             = realloc(rim->current_Ls, sizeof(int)*(N-1));
        rim->encounter_map          = realloc(rim->encounter_map,sizeof(int)*N);

        // Only need this stuff for Listing 3
        // rim->particles_backup_try   = realloc(rim->particles_backup_try,sizeof(struct reb_particle)*N);
        // rim->f0                     = realloc(rim->f0,sizeof(double)*((N-1)*(N-2))/2);
        // rim->f0_peris               = realloc(rim->f0_peris,sizeof(double)*(N-1));

        rim->allocatedN = N;

        //rim->close_encounters   = realloc(rim->close_encounters,sizeof(int)*N); // TLu
    }
    if (rim->safe_mode || rim->recalculate_coordinates_this_timestep){
        if (rim->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r);
            reb_warning(r,"MERCURIUS: Recalculating heliocentric coordinates but coordinates were not synchronized before.");
        }
        reb_integrator_mercurius_inertial_to_dh(r);
        rim->recalculate_coordinates_this_timestep = 0;
    }

    if (rim->recalculate_dcrit_this_timestep){
        rim->recalculate_dcrit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r);
            reb_integrator_mercurius_inertial_to_dh(r);
            rim->recalculate_coordinates_this_timestep = 0;
            reb_warning(r,"MERCURIUS: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        rim->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        for (int i=1;i<N;i++){
            rim->dcrit[i] = reb_integrator_mercurius_calculate_dcrit_for_particle(r, i);
        }
    }

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_MERCURIUS){
        reb_warning(r,"Mercurius has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_MERCURIUS;
    rim->mode = 0;

    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurius_L_mercury;
    }
}

void reb_integrator_mercurius_whfast_step(struct reb_simulation* const r, double dt){
  struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
  // If we are running WHFast K, L are always 0
  const int N =  r->N;
  for (int i = 1; i < r->N; i++){
    rim->current_Ls[i-1] = 0; // 0 means no close encounter
  }
  for (int i = 0; i < (N-1)*(N-2)/2; i++){
    rim->current_Ks[i] = 0; // 0 means no close encounter
  }

  reb_update_acceleration(r);
  if (rim->is_synchronized){
      reb_integrator_mercurius_interaction_step(r,dt/2.);
  }else{
      reb_integrator_mercurius_interaction_step(r,dt);
  }
  reb_integrator_mercurius_jump_step(r,dt/2.);
  reb_integrator_mercurius_com_step(r,dt);

  // WHFast Kepler step for all particles
  reb_integrator_mercurius_kepler_step(r, dt); // Kepler solver does NOT advance timestep

  reb_integrator_mercurius_jump_step(r,dt/2.);

  rim->is_synchronized = 0;
  if (rim->safe_mode){
      reb_integrator_mercurius_synchronize(r);
  }
}


void reb_integrator_mercurius_bs_step(struct reb_simulation* const r, double dt){
  struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
  reb_update_acceleration(r);
  if (rim->is_synchronized){
      reb_integrator_mercurius_interaction_step(r, dt/2.); // Pdot for B
  }else{
      reb_integrator_mercurius_interaction_step(r, dt);
  }
  reb_integrator_mercurius_jump_step(r,dt/2.);
  reb_integrator_mercurius_com_step(r,dt);

  reb_mercurius_encounter_predict(r);
  reb_mercurius_bs_encounter_step(r, dt);
  reb_integrator_mercurius_jump_step(r,dt/2.);

  rim->is_synchronized = 0;
  if (rim->safe_mode){
      reb_integrator_mercurius_synchronize(r); // Interaction here: PDot for B
  }
}

int pindex(int i, int j, int N){
  return (i-1)*N-((i-1)*(2+i)/2)+j-i-1;
}

// In listing 2 can use this for errything
int effF(struct reb_simulation* const r){
  struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
  const int N = r->N;
  const double* const dcrit = rim->dcrit;
  const double peri = rim->peri;
  int c0 = 0;

  // This is inefficient for many test particles...
  for (int i = 0; i < N; i++){
    for (int j = i + 1; j < N; j++){
      // Needed for both
      const double dx = r->particles[i].x - r->particles[j].x;
      const double dy = r->particles[i].y - r->particles[j].y;
      const double dz = r->particles[i].z - r->particles[j].z;
      const double d = sqrt(dx*dx + dy*dy + dz*dz);

      // Naive condition
      double dcritmax = MAX(dcrit[i],dcrit[j]);
      dcritmax *= 1.21;

      double fcond = d - dcritmax;
      double fcond_peri = d - peri;

/*
      // Velocity dependent condition

      const double vx = r->particles[i].vx - r->particles[j].vx;
      const double vy = r->particles[i].vy - r->particles[j].vy;
      const double vz = r->particles[i].vz - r->particles[j].vz;
      const double v = sqrt(vx*vx + vy*vy + vz*vz);

      double F_vel = d / (sqrt(3. * v * v + r->G * (r->particles[i].m + r->particles[j].m) / d));
      double fcond = F_vel - 30. * r->dt; // velocity dependent condition
      double fcond_peri = fcond;
*/


      // Body-body
      if (i != 0 && fcond < 0.0){
        rim->current_Ks[pindex(i,j,N)] = 1;
        c0 = 1;
      }

      // Check for close encounter with central body
      if (i == 0 && fcond_peri < 0.0){
        rim->current_Ls[j-1] = 1;
        c0 = 1;
      }
    }
  }
  return c0;
}

// Check F for pre-timestep using backup particles
int Fp(struct reb_simulation* const r, int* L_rej){
  struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
  const int N = r->N;
  const double* const dcrit = rim->dcrit;
  const double peri = rim->peri;
  int c0 = 0;

  for (int i = 0; i < N; i++){
    for (int j = i + 1; j < N; j++){
      // pre-timestep
      const double dx0 = rim->particles_backup[i].x - rim->particles_backup[j].x;
      const double dy0 = rim->particles_backup[i].y - rim->particles_backup[j].y;
      const double dz0 = rim->particles_backup[i].z - rim->particles_backup[j].z;
      const double d0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

      double dcritmax = MAX(dcrit[i],dcrit[j]);
      dcritmax *= 1.21;

      double f0 = d0 - dcritmax;
      double f0_peri = d0 - peri;

      // Body-body
      if (i != 0){
        // Populate f0 array so don't need to recalculate
        rim->f0[pindex(i,j,N)] = f0; // probably needs to be fact checked, but seems to be working for now
        if (f0 < 0.0){
          rim->current_Ks[pindex(i,j,N)] = 1;
          c0 = 1;
        }
      }

      // Check for close encounter with central body
      if (i == 0){
          rim->f0_peris[j-1] = f0_peri;
          if (f0_peri < 0.0){
            rim->current_Ls[j-1] = 1;
            c0 = 1;
            *L_rej = 1;
        }
      }
    }
  }

  return c0;
}

// Checks Ftry both pre- and post- timestep as in Listing 3 of Hernandez & Dehnen (2023)
// These do not alter K inherently
int F_cond(struct reb_simulation* const r, int cond){
  struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
  const int N = r->N;
  const double* const dcrit = rim->dcrit;
  const double peri = rim->peri;
  int ctry = 0;

  // Ftry only assesses planet-planet
  for (int i = 1; i < N; i++){
    for (int j = i + 1; j < N; j++){
      const double dx = r->particles[i].x - r->particles[j].x;
      const double dy = r->particles[i].y - r->particles[j].y;
      const double dz = r->particles[i].z - r->particles[j].z;
      const double d = sqrt(dx*dx + dy*dy + dz*dz);

      double f0 = rim->f0[pindex(i,j,N)];;
      double f = d - peri;;

      int c0;
      int c;
      if (cond == 0){ // this is Ftry
        c0 = f0 > 0;
        c = (f0 + f) > 0;

        if (c0 != c){
          ctry = 1;
        }
      }

      else{ // this is Falt
        c0 = f0 < 0;
        c = (f0 + f) < 0;

        if (c0 == c){
          ctry = 1;
        }
      }
    }
  }
  return ctry;
}

// This is Listing 2
void reb_integrator_mercurius_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    const int N = r->N;
    const double* const dcrit = rim->dcrit;

    // Make copy of particles
    memcpy(rim->particles_backup,r->particles,N*sizeof(struct reb_particle));

    // set current L at the beginning of each timestep
    for (int i = 1; i < N; i++){
      rim->current_Ls[i-1] = 0;
    }

    for (int i = 0; i < (N-1)*(N-2)/2; i++){
      rim->current_Ls[i-1] = 0;
    }

    int c0 = effF(r); // Check initial condition

    if (!c0){
      reb_integrator_mercurius_whfast_step(r, r->dt);

      int fy1 = effF(r);

      if (fy1){
        // This means close encounter post-timestep but NOT pre-timestep: we reject the WHFast step and integrate the backup particles with IAS15
        for (int i=0; i<N; i++){
            r->particles[i] = rim->particles_backup[i];
        }
        reb_integrator_mercurius_bs_step(r, r->dt);
      }
    }
    else{
      reb_integrator_mercurius_bs_step(r, r->dt);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

// TLu part2 from Hernandez & Dehnen 2023 Listing 3
// simple version: not particle-wise. Integrate whole sim w/ WHFast vs BS
/*
void reb_integrator_mercurius_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    const int N = r->N;
    const double* const dcrit = rim->dcrit;

    for (int i = 1; i < N; i++){
      rim->current_Ls[i-1] = 0; // 0 means no close encounter
    }
    rim->current_K = 0;

    // Make copy of particles
    memcpy(rim->particles_backup,r->particles,N*sizeof(struct reb_particle));

    // ytry for each particle
    int L_rej = 0;
    int c0 = Fp(r, &L_rej);

    // If there has been a close encounter, advance whole sim with BS
    if (c0){
      reb_integrator_mercurius_bs_step(r, r->dt);
    }
    // Otherwise, advance with WHFast
    else{
      reb_integrator_mercurius_whfast_step(r, r->dt);
    }

    // If close encounter with sun was detected, immediately accept BS step
    // otherwise:
    if (L_rej == 0){ // may be a better way to do this
      // Backup particles from first try in case alternative integrator is rejected
      memcpy(rim->particles_backup_try,r->particles,N*sizeof(struct reb_particle));

      // Now we assess if the step should be accepted by evaluating Ftry
      int ctry = F_cond(r,0);
      if (ctry){
        // Reject first integrator, try alternative integrators
        // In both cases, reset to backup values
        for (int i=0; i<N; i++){
            r->particles[i] = rim->particles_backup[i];
        }

        if (c0){
          // Means we tried BS and rejected. Now try with WHFast
          reb_integrator_mercurius_whfast_step(r, r->dt);
        }
        else{
          rim->current_K = 1;
          reb_integrator_mercurius_bs_step(r, r->dt);
        }
        // Now assess Falt
        int calt = F_cond(r,1);
          if (calt){
            // Reject alternative integrator, use original integrators
            for (int i=0; i<N; i++){
                r->particles[i] = rim->particles_backup_try[i];
            }
          }
        }
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}
*/

void reb_integrator_mercurius_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    if (rim->is_synchronized == 0){
        r->gravity = REB_GRAVITY_MERCURIUS; // needed here again for SimulationArchive
        rim->mode = 0;
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurius_L_mercury;
        }
        //printf("Before interaction 2\n");
        reb_update_acceleration(r);
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);

        reb_integrator_mercurius_dh_to_inertial(r);

        rim->recalculate_coordinates_this_timestep = 1;
        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurius_reset(struct reb_simulation* r){
    r->ri_mercurius.L = NULL;
    r->ri_mercurius.mode = 0;
    r->ri_mercurius.encounterN = 0;
    r->ri_mercurius.encounterNactive = 0;
    r->ri_mercurius.hillfac = 4; // TLu changed to Hernandez (2023)
    //r->ri_mercurius.peri = 0.; // TLu changed to Hernandez (2023)
    r->ri_mercurius.tponly_encounter = 0;
    r->ri_mercurius.recalculate_coordinates_this_timestep = 0;
    // Internal arrays (only used within one timestep)
    free(r->ri_mercurius.particles_backup);
    r->ri_mercurius.particles_backup = NULL;
    free(r->ri_mercurius.particles_backup_additionalforces);
    r->ri_mercurius.particles_backup_additionalforces = NULL;
    free(r->ri_mercurius.encounter_map);
    r->ri_mercurius.encounter_map = NULL;
    r->ri_mercurius.allocatedN = 0;
    r->ri_mercurius.allocatedN_additionalforces = 0;
    // dcrit array
    free(r->ri_mercurius.dcrit);
    r->ri_mercurius.dcrit = NULL;
    r->ri_mercurius.dcrit_allocatedN = 0;

    free(r->ri_mercurius.particles_backup_try);
    r->ri_mercurius.particles_backup_try = NULL;

    free(r->ri_mercurius.f0);
    r->ri_mercurius.f0 = NULL;

    free(r->ri_mercurius.f0_peris);
    r->ri_mercurius.f0_peris = NULL;

    //free(r->ri_mercurius.current_Ls);
    //r->ri_mercurius.current_Ls = NULL;

    // free(r->ri_mercurius.close_encounters);
    // r->close_encounters = NULL;
}
/*
void reb_integrator_tlu_bs_step(struct reb_simulation* r, double dt){
  int nRefineMax = 8;
  double T[0][nRefineMax][nRefineMax] ={0};
  double T[0][nRefineMax] ={0};
}
*/
