/**
 * @file    integrator_trace.c
 * @brief   trace, a modified version of John Chambers' MERCURY algorithm
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
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_trace.h"
#include "integrator_whfast.h"
#include "integrator_bs.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

double reb_integrator_trace_switch_default(struct reb_simulation* const r, int i, int j){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  double dcriti = 0.0;
  double dcritj = 0.0;

  const double m0 = r->particles[0].m;

  if (r->particles[i].m != 0){

    const double dxi  = r->particles[i].x;  // in dh
    const double dyi  = r->particles[i].y;
    const double dzi  = r->particles[i].z;
    const double dvxi = r->particles[i].vx - r->particles[0].vx;
    const double dvyi = r->particles[i].vy - r->particles[0].vy;
    const double dvzi = r->particles[i].vz - r->particles[0].vz;
    const double _ri = sqrt(dxi*dxi + dyi*dyi + dzi*dzi);
    const double v2i = dvxi*dvxi + dvyi*dvyi + dvzi*dvzi;

    const double GMi = r->G*(m0+r->particles[i].m);
    const double ai = GMi*_ri / (2.*GMi - _ri*v2i);

    //struct reb_orbit o1 = reb_tools_particle_to_orbit(r->G, r->particles[i], r->particles[0]);
    //const double ai = o1.a;
    dcriti = ri_tr->hillfac*ai*cbrt(r->particles[i].m/(3.*m0));
  }

  if (r->particles[j].m != 0){

    const double dxj  = r->particles[j].x;  // in dh
    const double dyj  = r->particles[j].y;
    const double dzj  = r->particles[j].z;
    const double dvxj = r->particles[j].vx - r->particles[0].vx;
    const double dvyj = r->particles[j].vy - r->particles[0].vy;
    const double dvzj = r->particles[j].vz - r->particles[0].vz;
    const double _rj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);
    const double v2j = dvxj*dvxj + dvyj*dvyj + dvzj*dvzj;

    const double GMj = r->G*(m0+r->particles[j].m);
    const double aj = GMj*_rj / (2.*GMj - _rj*v2j);

    //struct reb_orbit o2 = reb_tools_particle_to_orbit(r->G, r->particles[j], r->particles[0]);
    //const double aj = o2.a;

    dcritj = ri_tr->hillfac*aj*cbrt(r->particles[j].m/(3.*m0));
  }

  const double dx = r->particles[i].x - r->particles[j].x;
  const double dy = r->particles[i].y - r->particles[j].y;
  const double dz = r->particles[i].z - r->particles[j].z;
  const double d = sqrt(dx*dx + dy*dy + dz*dz);

  // This stuff seems unimportant...
  /*
  // Interaction hamiltonian
  const double hi = (r->G * r->particles[j].m * r->particles[i].m) / d;

  // Kepler Hamiltonian i
  const double hki = (v2i / (r->particles[i].m * r->particles[i].m)) / (2 * r->particles[i].m) - ((r->G * r->particles[0].m * r->particles[i].m) / _ri);
  const double hkj = (v2j / (r->particles[j].m * r->particles[j].m)) / (2 * r->particles[j].m) - ((r->G * r->particles[0].m * r->particles[j].m) / _ri);
  //if (ri_tr->print){
  //  printf("Hamiltonians at %f for %d %d: %e %e\n", r->t, i, j, hi/hki, hi/hkj);
  //}

  if ((ai <= 0.0 || aj <= 0.0) && (hi/hki > 5e-15 || hi/hkj > 5e-15)){
    if (hi/hki > 5e-15 || hi/hkj > 5e-15){
      return -1.0;
    }
    else{
      return 1.0;
    }
  }
  */

  // Use traditional switching function
  double dcritmax = MAX(dcriti,dcritj);
  dcritmax *= 1.21;

  double fcond = d - dcritmax;
  return fcond;
}

double reb_integrator_trace_peri_switch_default(const struct reb_simulation* const r, int j){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const double peri = ri_tr->peri;

  const double dx = r->particles[0].x - r->particles[j].x;
  const double dy = r->particles[0].y - r->particles[j].y;
  const double dz = r->particles[0].z - r->particles[j].z;
  const double d = sqrt(dx*dx + dy*dy + dz*dz);

  double fcond_peri = d - peri;
  return fcond_peri;
}

double reb_integrator_trace_switch_fdot_peri(const struct reb_simulation* const r, int j){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const double vfacp = ri_tr->vfac_p;

  const double dx = r->particles[j].x;
  const double dy = r->particles[j].y;
  const double dz = r->particles[j].z;
  const double d2 = dx*dx + dy*dy + dz*dz;
  const double d = sqrt(d2);
  //printf("%d %f\n", j, d);

  const double dvx = r->particles[j].vx - r->particles[0].vx;
  const double dvy = r->particles[j].vy - r->particles[0].vy;
  const double dvz = r->particles[j].vz - r->particles[0].vz;
  const double v2 = dvx * dvx + dvy * dvy + dvz * dvz;

  const double hx = (dy*dvz - dz*dvy); 					// specific angular momentum vector
  const double hy = (dz*dvx - dx*dvz);
  const double hz = (dx*dvy - dy*dvx);
  const double h = sqrt ( hx*hx + hy*hy + hz*hz );

  // This only works for bound orbits!
  const double fdot = h / (d2);
  const double peff = (2 * M_PI / fdot); // effective period
  double fcond_peri = (peff / r->dt) - vfacp;

  // Failsafe - use velocity dependent condition
  double f_vel = d / sqrt(3. * v2 + (r->G * (r->particles[0].m + r->particles[j].m) / d));
  double fcond_vel = (f_vel/r->dt) - 1.;

  // if one is violated both are
  double fcond = MIN(fcond_peri, fcond_vel);
  //if (fcond < 0.0){
  //  printf("%f %d %f %f\n", r->t, j, fcond_peri, fcond_vel);
  //}
  return fcond;
}

double reb_integrator_trace_switch_vdiff_peri(const struct reb_simulation* const r, int j){
  struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[j], r->particles[0]);
  if (o.a < 0.0){
    return 1.; //unbound orbits have no pericenter approach
  }
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const double vfacp = ri_tr->vfac_p;

  const double dvx = r->particles[j].vx;
  const double dvy = r->particles[j].vy;
  const double dvz = r->particles[j].vz;
  const double dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

  const double vcirc = sqrt(r->G * (r->particles[0].m + r->particles[j].m) / o.a);
  const double vdiff = dv / vcirc;
  double fcond_peri = vfacp - vdiff;

  return fcond_peri;
}

void reb_integrator_trace_inertial_to_dh(struct reb_simulation* r){
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
    r->ri_tr.com_pos = com_pos;
    r->ri_tr.com_vel = com_vel;
}

void reb_integrator_trace_dh_to_inertial(struct reb_simulation* r){
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
    particles[0].x = r->ri_tr.com_pos.x - temp.x;
    particles[0].y = r->ri_tr.com_pos.y - temp.y;
    particles[0].z = r->ri_tr.com_pos.z - temp.z;

    for (int i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += r->ri_tr.com_vel.x;
        particles[i].vy += r->ri_tr.com_vel.y;
        particles[i].vz += r->ri_tr.com_vel.z;
    }
    particles[0].vx = r->ri_tr.com_vel.x - temp.vx;
    particles[0].vy = r->ri_tr.com_vel.y - temp.vy;
    particles[0].vz = r->ri_tr.com_vel.z - temp.vz;
}

void reb_integrator_trace_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    reb_update_acceleration(r);
    for (int i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_trace_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;

    struct reb_simulation_integrator_trace* ri_tr = &(r->ri_tr);
    const int current_L = ri_tr->current_L;

    const int N_active = r->N_active==-1?r->N:r->N_active;

    // If TP type 1, use r->N. Else, use N_active.
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
    for (int i=1;i<N_all;i++){
        particles[i].x += dt*px*(1-current_L);
        particles[i].y += dt*py*(1-current_L);
        particles[i].z += dt*pz*(1-current_L);
    }
}

void reb_integrator_trace_com_step(struct reb_simulation* const r, double dt){
    r->ri_tr.com_pos.x += dt*r->ri_tr.com_vel.x;
    r->ri_tr.com_pos.y += dt*r->ri_tr.com_vel.y;
    r->ri_tr.com_pos.z += dt*r->ri_tr.com_vel.z;
}

// Old Kepler
void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    //struct reb_particle* restrict const particles = r->particles;
    struct reb_simulation_integrator_trace* ri_tr = &(r->ri_tr);
    const int N = r->N;
    for (int i=1;i<N;i++){
        //if (ri_tr->encounter_map[i] != 0){
          reb_whfast_kepler_solver(r,r->particles,r->G*r->particles[0].m,i,dt); // in dh
        //}
    }
}

void reb_integrator_trace_bs_step(struct reb_simulation* const r, const double _dt){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);

  if (ri_tr->encounterN < 2){
    // No close encounters, skip
    return;
  }

  int i_enc = 0;
  ri_tr->encounterNactive = 0;
  for (unsigned int i=0; i<r->N; i++){
      if(ri_tr->encounter_map_internal[i]){
          struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
          r->particles[i] = ri_tr->particles_backup_try[i]; // Coordinates before WHFast step, overwrite particles with close encounters
          ri_tr->encounter_map[i_enc] = i;
          i_enc++;
          if (r->N_active==-1 || i<r->N_active){
              ri_tr->encounterNactive++;
              if (ri_tr->tponly_encounter){
                    ri_tr->particles_backup_try[i] = tmp;         // Make copy of particles after the kepler step.
                                                            // used to restore the massive objects' states in the case
                                                            // of only massless test-particle encounters
                }
          }
      }
  }
  // printf("Encounter step triggered %f\n", r->t);
  //exit(1);

  ri_tr->mode = 1;
  // run
  const double old_dt = r->dt;
  const double old_t = r->t;
  double t_needed = r->t + _dt;
  //reb_integrator_bs_reset(r);

  r->dt = _dt; // start with a small timestep.
  //printf("BEGIN!!! %f %f ", r->t, r->dt);

  while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 ){
    // In case of overshoot
    if (r->t+r->dt >  t_needed){
      r->dt = t_needed-r->t;
    }

    struct reb_particle star = r->particles[0]; // backup velocity
    r->particles[0].vx = 0; // star does not move in dh
    r->particles[0].vy = 0;
    r->particles[0].vz = 0;

    reb_integrator_bs_part2(r);

    const double dx = r->particles[2].x;
    const double dy = r->particles[2].y;
    const double dz = r->particles[2].z;
    const double d = sqrt(dx*dx + dy*dy + dz*dz);

    reb_collision_search(r);

    // Now, r->dt is the proposed next step
    r->particles[0].vx = star.vx; // restore every timestep for collisions
    r->particles[0].vy = star.vy;
    r->particles[0].vz = star.vz;

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

  // if only test particles encountered massive bodies, reset the
  // massive body coordinates to their post Kepler step state
  if(ri_tr->tponly_encounter){
      for (unsigned int i=1; i < ri_tr->encounterNactive; i++){
          unsigned int mi = ri_tr->encounter_map[i];
          r->particles[mi] = ri_tr->particles_backup_try[mi];
      }
  }

  //for (unsigned int i = 0; i<r->N; i++){
  //  printf("%f %d %e %e %e %e %e\n", r->t, i, r->particles[i].m, r->particles[i].x, r->particles[i].y, r->particles[i].vx, r->particles[i].vy);
  //}

  r->t = old_t;
  r->dt = old_dt;
  ri_tr->mode = 0;
  // return reject;
}

void reb_integrator_trace_kepler_step(struct reb_simulation* const r, const double _dt){
  int rej = 0;
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  memcpy(ri_tr->particles_backup_try,r->particles,r->N*sizeof(struct reb_particle));
  reb_integrator_trace_whfast_step(r, _dt);
  reb_integrator_trace_bs_step(r, _dt);
  // return rej;
}


void reb_integrator_trace_part1(struct reb_simulation* r){
   //printf("TRACE part 1\n");
    if (r->var_config_N){
        reb_warning(r,"TRACE does not work with variational equations.");
    }

    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    const int N = r->N;

    if (ri_tr->allocatedN<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_tr->particles_backup       = realloc(ri_tr->particles_backup,sizeof(struct reb_particle)*N);
        ri_tr->current_Ks             = realloc(ri_tr->current_Ks, sizeof(int*)*N); // This is inefficient for now, can be Nactive instead of N
        for (int k = 0; k < N; ++k) {
            ri_tr->current_Ks[k]      = realloc(ri_tr->current_Ks[k], sizeof(int)*N);
        }

        ri_tr->encounter_map          = realloc(ri_tr->encounter_map,sizeof(int)*N);
        ri_tr->encounter_map_internal = realloc(ri_tr->encounter_map_internal,sizeof(int)*N); // Do we need this now?

        // Only need this stuff for Listing 3
        ri_tr->particles_backup_try   = realloc(ri_tr->particles_backup_try,sizeof(struct reb_particle)*N);
        ri_tr->allocatedN = N;
    }

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"TRACE only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_TRACE){
        reb_warning(r,"TRACE has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }

    // Don't need to do this every timestep... fix this later
    reb_integrator_trace_inertial_to_dh(r);

    // Switching functions
    if (ri_tr->S == NULL){
      ri_tr->S = reb_integrator_trace_switch_default;
    }

    if (ri_tr->S_peri == NULL){
      ri_tr->S_peri = reb_integrator_trace_switch_fdot_peri;
    }

    r->gravity = REB_GRAVITY_TRACE;
    ri_tr->mode = 0;
    ri_tr->collision = 0;

    // Clear encounter maps
    for (unsigned int i=0; i<r->N; i++){
      ri_tr->encounter_map[i] = 0;
      ri_tr->encounter_map_internal[i] = 0;
    }
    ri_tr->encounter_map_internal[0] = 1;
}

void reb_integrator_trace_F_start(struct reb_simulation* const r){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const int Nactive = r->N_active==-1?r->N:r->N_active;

  if (r->testparticle_type == 1){
      ri_tr->tponly_encounter = 0; // testparticles affect massive particles
  }else{
      ri_tr->tponly_encounter = 1;
  }

  // Switching functions
  double (*_switch) (const struct reb_simulation* const r, unsigned int i, unsigned int j) = r->ri_tr.S;
  double (*_switch_peri) (const struct reb_simulation* const r, unsigned int j) = r->ri_tr.S_peri;

  // Check for pericenter CE
  // test particles cannot have pericenter CEs
  for (int j = 1; j < Nactive; j++){
    double fcond_peri = _switch_peri(r, j);
    if (fcond_peri < 0.0){
      ri_tr->current_L = 1;
      //if (ri_tr->print){
      //printf("Flagged %d peri approach at %f %f\n", j, r->t, fcond_peri);
      //}
      if (j < Nactive){ // Two massive particles have a close encounter
          ri_tr->tponly_encounter = 0;
      }
    }
  }

  // Body-body
  // there cannot be TP-TP CEs
  for (int i = 1; i < Nactive; i++){
    for (int j = i + 1; j < N; j++){

      double fcond = _switch(r, i, j);

      if (fcond < 0.0){
        //if (ri_tr->print){
        //printf("Flagged %d %d CE at %f %f\n", i, j, r->t, fcond);
        //}

        if (ri_tr->encounter_map_internal[i] == 0){
            ri_tr->encounter_map_internal[i] = i;
            ri_tr->encounterN++;
        }
        if (ri_tr->encounter_map_internal[j] == 0){
            ri_tr->encounter_map_internal[j] = j;
            ri_tr->encounterN++;
        }

        if (j < Nactive){ // Two massive particles have a close encounter
            ri_tr->tponly_encounter = 0;
        }

        ri_tr->current_Ks[i][j] = 1;
      }
    }
  }
}

int reb_integrator_trace_Fcond(struct reb_simulation* const r){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const int Nactive = r->N_active==-1?r->N:r->N_active;
  //int print = ri_tr->print;

  int new_c = 0; // New CEs

  // Switching functions
  double (*_switch) (const struct reb_simulation* const r, unsigned int i, unsigned int j) = r->ri_tr.S;
  double (*_switch_peri) (const struct reb_simulation* const r, unsigned int j) = r->ri_tr.S_peri;

  if (r->testparticle_type == 1){
      ri_tr->tponly_encounter = 0; // testparticles affect massive particles
  }else{
      ri_tr->tponly_encounter = 1;
  }

  // Check for pericenter CE
  // test particles cannot have pericenter CEs

  for (int j = 1; j < Nactive; j++){
    double fcond_peri = _switch_peri(r, j);
    if (fcond_peri < 0.0 && ri_tr->current_L == 0){
      ri_tr->current_L = 1;
      new_c = 1;

      if (j < Nactive){ // Two massive particles have a close encounter
          ri_tr->tponly_encounter = 0;
      }
    }
  }

  // Body-body
  // there cannot be TP-TP CEs
  for (int i = 1; i < Nactive; i++){
    for (int j = i + 1; j < N; j++){

      double fcond = _switch(r, i, j);

      if (fcond < 0.0){
        if (ri_tr->encounter_map_internal[i] == 0){
            ri_tr->encounter_map_internal[i] = i;
            ri_tr->encounterN++;
        }
        if (ri_tr->encounter_map_internal[j] == 0){
            ri_tr->encounter_map_internal[j] = j;
            ri_tr->encounterN++;
        }

        if (j < Nactive){ // Two massive particles have a close encounter
            ri_tr->tponly_encounter = 0;
        }

        // Checks for switching Kij 0->1. Initialized as all 0 the first time of asking.
        if (ri_tr->current_Ks[i][j] == 0){
          ri_tr->current_Ks[i][j] = 1;
          new_c = 1;
          //if (ri_tr->print){
          //  printf("CE at %f detected between %d %d\n", r->t,i, j);
          //}
        }
      }
    }
  }

  return new_c;
}

int reb_integrator_trace_F_check(struct reb_simulation* const r, int old_N){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const int Nactive = r->N_active==-1?r->N:r->N_active;
  int reject = 0;

  // Switching functions
  double (*_switch) (const struct reb_simulation* const r, unsigned int i, unsigned int j) = r->ri_tr.S;
  double (*_switch_peri) (const struct reb_simulation* const r, unsigned int j) = r->ri_tr.S_peri;

  // Check for L 0 -> 1
  // Deal with collisions here later
  if (ri_tr->current_L == 0){
    for (int j = 1; j < Nactive; j++){
      double fcond_peri = _switch_peri(r, j);
      if (fcond_peri < 0.0){
        ri_tr->current_L = 1;
        reject = 1;

        //printf("REJECT BC %d peri approach at %f %f\n", j, r->t, fcond_peri);
        break;
      }
    }
  }


  // Check for K 0 -> 1
  for (int i = 1; i < Nactive; i++){
    for (int j = i + 1; j < N; j++){
      if (ri_tr->current_Ks[i][j] == 0){
        double fcond = _switch(r, i, j);
        if (fcond < 0.0){
          // Need to fix internal maps...
          if (ri_tr->encounter_map_internal[i] == 0){
              ri_tr->encounter_map_internal[i] = i;
              ri_tr->encounterN++;
          }
          if (ri_tr->encounter_map_internal[j] == 0){
              ri_tr->encounter_map_internal[j] = j;
              ri_tr->encounterN++;
          }
          if (j < Nactive){ // Two massive particles have a close encounter
              ri_tr->tponly_encounter = 0;
          }
          ri_tr->current_Ks[i][j] = 1;
          reject = 1;
        }
      }
    }
  }
  return reject;
}

// This is Listing 2
void reb_integrator_trace_part2(struct reb_simulation* const r){
  //printf("TRACE part 2\n");
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    const int N = r->N;
    //if (r->t >= 6.54e5 * 2. * M_PI){
      //ri_tr->print = 1;
    //}
    // Make copy of particles
    memcpy(ri_tr->particles_backup,r->particles,N*sizeof(struct reb_particle));
    ri_tr->encounterN = 1;
    ri_tr->current_L = 0;

    //printf("Initial: %f\n", r->particles[2].vx);

    for (int i = 0; i < N; i++){
      for (unsigned int j = i + 1; j < N; j++){
          ri_tr->current_Ks[i][j] = 0;
      }
    }

    //reb_integrator_trace_F_start(r);
    int rej = reb_integrator_trace_Fcond(r); // output means nothing here
    if (ri_tr->current_L){ //more efficient way to check if we need to redo this...
      // Pericenter close encounter detected. We integrate the entire simulation with BS
      //printf("Flagged Peri %f\n", r->t);
      ri_tr->encounter_map_internal[0] = 1;
      ri_tr->encounterN = N;
      for (int i = 1; i < N; i++){
        ri_tr->encounter_map_internal[i] = i; // Identity map
      }
      ri_tr->encounterNactive = ((r->N_active==-1)?r->N:r->N_active);
    }


    reb_integrator_trace_interaction_step(r, r->dt/2.);
    reb_integrator_trace_jump_step(r, r->dt/2.);
    reb_integrator_trace_kepler_step(r, r->dt); // always accept this
    reb_integrator_trace_com_step(r,r->dt);
    reb_integrator_trace_jump_step(r, r->dt/2.);
    reb_integrator_trace_interaction_step(r, r->dt/2.);

    // Check for new close_encounters
    //if (ri_tr->print){
    //  printf("\nREJECTION CHECK\n");
    //}
    //if (reb_integrator_trace_F_check(r, N) && !ri_tr->collision){
    if (reb_integrator_trace_Fcond(r)){
      // REJECT STEP
      // reset simulation and try again with new timestep
      for (int i=0; i<N; i++){
          // Reject & reset
          r->particles[i] = ri_tr->particles_backup[i];
      }

      if (ri_tr->current_L){ //more efficient way to check if we need to redo this...
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        //printf("REJECT Flagged Peri %f\n", r->t);
        ri_tr->encounter_map_internal[0] = 1;
        ri_tr->encounterN = N;
        for (int i = 1; i < N; i++){
          ri_tr->encounter_map_internal[i] = i; // Identity map
        }
        ri_tr->encounterNactive = ((r->N_active==-1)?r->N:r->N_active);
      }


      reb_integrator_trace_interaction_step(r, r->dt/2.);
      reb_integrator_trace_jump_step(r, r->dt/2.);
      reb_integrator_trace_kepler_step(r, r->dt); // always accept this
      reb_integrator_trace_com_step(r,r->dt);
      reb_integrator_trace_jump_step(r, r->dt/2.);
      reb_integrator_trace_interaction_step(r, r->dt/2.);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
    ri_tr->mode = 0;
    r->gravity = REB_GRAVITY_TRACE; // Is this needed?
    reb_integrator_trace_dh_to_inertial(r);
}

void reb_integrator_trace_synchronize(struct reb_simulation* r){
}

void reb_integrator_trace_reset(struct reb_simulation* r){
    r->ri_tr.mode = 0;
    r->ri_tr.encounterN = 0;
    r->ri_tr.encounterNactive = 0;
    r->ri_tr.hillfac = 4; // TLu changed to Hernandez (2023)
    r->ri_tr.vfac_p = 32.;

    //r->ri_tr.peri = 0.; // TLu changed to Hernandez (2023)
    // Internal arrays (only used within one timestep)
    free(r->ri_tr.particles_backup);
    r->ri_tr.particles_backup = NULL;

    free(r->ri_tr.encounter_map);
    r->ri_tr.encounter_map = NULL;
    r->ri_tr.allocatedN = 0;
    r->ri_tr.allocatedN_additionalforces = 0;

    free(r->ri_tr.particles_backup_try);
    r->ri_tr.particles_backup_try = NULL;

    if (r->ri_tr.current_Ks){
        for (int k = 0; k < r->N; ++k) {
            r->ri_tr.current_Ks[k] = NULL;
        }
        free(r->ri_tr.current_Ks);
        r->ri_tr.current_Ks = NULL;
    }

    r->ri_tr.current_L = 0;

    free(r->ri_tr.encounter_map_internal);
    r->ri_tr.encounter_map_internal = NULL;

    r->ri_tr.S = NULL;
    r->ri_tr.S_peri = NULL;

}
