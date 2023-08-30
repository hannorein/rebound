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
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_trace.h"
#include "integrator_whfast.h"
#include "integrator_mercurius.h"
#include "integrator_bs.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

// Can we just use Mercurius for these? Maybe not if Jacobi...

struct reb_vec3d reb_integrator_trace_switch_default(struct reb_simulation* const r, int i, int j){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  struct reb_vec3d a = {0};
  // const double* const dcrit = ri_tr->dcrit;
  const double peri = ri_tr->peri;

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

    const double GM = r->G*(m0+r->particles[i].m);
    const double ai = GM*_ri / (2.*GM - _ri*v2i);
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

    const double GM = r->G*(m0+r->particles[j].m);
    const double aj = GM*_rj / (2.*GM - _rj*v2j);
    dcritj = ri_tr->hillfac*aj*cbrt(r->particles[j].m/(3.*m0));
  }

  const double dx = r->particles[i].x - r->particles[j].x;
  const double dy = r->particles[i].y - r->particles[j].y;
  const double dz = r->particles[i].z - r->particles[j].z;
  const double d = sqrt(dx*dx + dy*dy + dz*dz);

  // Use traditional switching function
  double dcritmax = MAX(dcriti,dcritj);
  dcritmax *= 1.21;

  double fcond = d - dcritmax;
  double fcond_peri = d - peri;

  a.x = fcond;
  a.y = fcond_peri;

  return a;
}

struct reb_vec3d reb_integrator_trace_switch_velocity(const struct reb_simulation* const r, int i, int j){
  struct reb_vec3d conds = {0};
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const double vfac = ri_tr->vfac;
  const double vfacp = ri_tr->vfac_p;

  const double dx = r->particles[i].x - r->particles[j].x;
  const double dy = r->particles[i].y - r->particles[j].y;
  const double dz = r->particles[i].z - r->particles[j].z;
  const double d = sqrt(dx*dx + dy*dy + dz*dz);

  const double vx = r->particles[i].vx - r->particles[j].vx;
  const double vy = r->particles[i].vy - r->particles[j].vy;
  const double vz = r->particles[i].vz - r->particles[j].vz;
  const double vsquared = vx*vx + vy*vy + vz*vz;

  const double F_vel = d / (sqrt(3. * vsquared + r->G * (r->particles[i].m + r->particles[j].m) / d));
  const double fcond = F_vel - vfac * r->dt;

  const double mu = r->G*(r->particles[i].m+r->particles[j].m);
  const double vcircsquared = mu/d;
  const double a = -mu/(vsquared - 2.*vcircsquared);

  const double vdiffsquared = vsquared - vcircsquared;
  const double vr = (dx*vx + dy*vy + dz*vz)/d;
  const double rvr = d*vr;
  const double muinv = 1./mu;
  const double ex = muinv*(vdiffsquared*dx - rvr*vx );
  const double ey = muinv*(vdiffsquared*dy - rvr*vy );
  const double ez = muinv*(vdiffsquared*dz - rvr*vz );
  const double e = sqrt( ex*ex + ey*ey + ez*ez );
  const double ome = 1. - e;

  // The fdot condition
  const double F_vel_peri = 2 * M_PI * sqrt(((ome*ome*ome)/(1+e))*(a*a*a/mu));
  double fcond_peri = F_vel_peri - vfacp * r->dt;

  conds.x = fcond;
  conds.y = fcond_peri;

  return conds;
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
        //if (i == 2){
        //printf("%e %e %e %e\n", r->t, particles[i].ax, particles[i].ay, particles[i].az);
      //}
    }
}

void reb_integrator_trace_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;

    struct reb_simulation_integrator_trace* ri_tr = &(r->ri_tr);
    const int current_L = ri_tr->current_L;

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
void reb_integrator_trace_kepler_step(struct reb_simulation* const r, double dt){
    //struct reb_particle* restrict const particles = r->particles;
    struct reb_simulation_integrator_trace* ri_tr = &(r->ri_tr);
    const int N = r->N;
    for (int i=1;i<N;i++){
        //if (ri_tr->encounter_map[i] != 0){
          reb_whfast_kepler_solver(r,r->particles,r->G*r->particles[0].m,i,dt); // in dh
        //}
    }
}

int reb_trace_bs_step(struct reb_simulation* const r, const double _dt){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);

  if (ri_tr->encounterN < 2){
    // No close encounters, skip
    return 0;
  }

  int reject = 0;

  int i_enc = 0;
  ri_tr->encounterNactive = 0;
  for (unsigned int i=0; i<r->N; i++){
      if(ri_tr->encounter_map[i]){
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

  ri_tr->mode = 1;
  // run
  const double old_dt = r->dt;
  const double old_t = r->t;
  double t_needed = r->t + _dt;
  //reb_integrator_bs_reset(r);

  r->dt = _dt; // start with a small timestep.
  //printf("BEGIN!!! %f %f ", r->t, r->dt);
  while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 ){
      struct reb_particle star = r->particles[0]; // backup velocity
      r->particles[0].vx = 0; // star does not move in dh
      r->particles[0].vy = 0;
      r->particles[0].vz = 0;
      //reb_update_acceleration(r);
      reb_integrator_bs_part2(r);
      //printf("TRACE proposed: %f %f %f\n", r->dt, r->t, t_needed);

      r->particles[0].vx = star.vx; // restore every timestep for collisions
      r->particles[0].vy = star.vy;
      r->particles[0].vz = star.vz;

      // TLUTS
      ri_tr->dt_proposed = r->dt;

      // Reject timestep
      if (r->dt < ri_tr->dt_shells[ri_tr->current_shell] && ri_tr->current_shell != 0){
        ri_tr->current_shell -= 1;
        reject = 1;
        //printf("Rejected, Reduce Timestep to %f ", ri_tr->dt_shells[ri_tr->current_shell]);
        break;
      }

      if (r->t+r->dt >  t_needed){
          r->dt = t_needed-r->t;
      }


      // Search and resolve collisions. Ignore for now

      reb_collision_search(r);

      // Do any additional post_timestep_modifications.
      // Note: post_timestep_modifications is called here but also
      // at the end of the full timestep. The function thus needs
      // to be implemented with care as not to do the same
      // modification multiple times. To do that, check the value of
      // r->ri_tr.mode
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


//    if(ri_tr->tponly_encounter){
//        for (int i=1;i<ri_tr->encounterNactive;i++){
//            unsigned int mi = ri_tr->encounter_map[i];
//            r->particles[mi] = ri_tr->particles_backup_try[mi];
//        }
//    }
    //printf("Fin\n");
    r->t = old_t;
    r->dt = old_dt;
    ri_tr->mode = 0;
    return reject;
    //exit(1);
}


void reb_integrator_trace_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"TRACE does not work with variational equations.");
    }

    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    const int N = r->N;
    //ri_tr->encounterN = N; // TLu temp for full IAS15 integration

    if (ri_tr->dcrit_allocatedN<N){
        // Need to safe these arrays in SimulationArchive
        //ri_tr->dcrit              = realloc(ri_tr->dcrit, sizeof(double)*N);
        ri_tr->dcrit_allocatedN = N; // better way around this maybe
        // If particle number increased (or this is the first step), need to calculate critical radii
        //ri_tr->recalculate_dcrit_this_timestep        = 1;
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
         ri_tr->recalculate_coordinates_this_timestep = 1;
    }
    if (ri_tr->allocatedN<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_tr->particles_backup       = realloc(ri_tr->particles_backup,sizeof(struct reb_particle)*N);

        ri_tr->current_Ks = realloc(ri_tr->current_Ks, sizeof(int)*((N-1)*(N-2))/2);
        ri_tr->delta_Ks   = realloc(ri_tr->delta_Ks, sizeof(int)*100);//sizeof(int)*((N-1)*(N-2))/2);
        //ri_tr->current_Ls             = realloc(ri_tr->current_Ls, sizeof(int)*(N-1));
        ri_tr->encounter_map          = realloc(ri_tr->encounter_map,sizeof(int)*N);
        ri_tr->encounter_map_backup   = realloc(ri_tr->encounter_map_backup,sizeof(int)*N);

        // Only need this stuff for Listing 3
        ri_tr->particles_backup_try   = realloc(ri_tr->particles_backup_try,sizeof(struct reb_particle)*N);
        ri_tr->allocatedN = N;

        //ri_tr->close_encounters   = realloc(ri_tr->close_encounters,sizeof(int)*N); // TLu
    }
    if (ri_tr->safe_mode || ri_tr->recalculate_coordinates_this_timestep){
        if (ri_tr->is_synchronized==0){
            reb_integrator_trace_synchronize(r);
            reb_warning(r,"TRACE: Recalculating heliocentric coordinates but coordinates were not synchronized before.");
        }
        reb_integrator_trace_inertial_to_dh(r);
        ri_tr->recalculate_coordinates_this_timestep = 0;
    }

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"trace only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_TRACE){
        reb_warning(r,"TRACE has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }

    if (ri_tr->S == NULL){
      ri_tr->S = reb_integrator_trace_switch_default;
    }

    if (r->t == 0){
      ri_tr->dt_proposed = r->dt;
    }

    if (ri_tr->dt_shells == NULL){
        // Allocate dt_shells
        ri_tr->dt_shells = realloc(ri_tr->dt_shells, sizeof(double)*ri_tr->nshells);
        for (int i = 0; i < ri_tr->nshells; i++){
          ri_tr->dt_shells[i] = pow(2., i) * r->dt;
        }
    }

    r->gravity = REB_GRAVITY_TRACE;
    ri_tr->mode = 0;
    ri_tr->print = 0;
}

// Particle-particle collision tracking. Explanation is in my notes.
int pindex(int i, int j, int N){
  return (i-1)*N-((i-1)*(2+i)/2)+j-i-1;
}

int F0(struct reb_simulation* const r){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const int Nactive = r->N_active==-1?r->N:r->N_active;
  struct reb_vec3d (*_switch) (const struct reb_simulation* const r, int i, int j) = r->ri_tr.S;
  int c0 = 0;

  for (int i = 0; i < Nactive; i++){
    for (int j = i + 1; j < N; j++){

      struct reb_vec3d cond = _switch(r, i, j);
      double fcond = cond.x;
      double fcond_peri = cond.y;

      // Check for close encounter with central body. No test particles
      if (i == 0 && j < Nactive && ri_tr->current_L == 0 && fcond_peri < 0.0){
        ri_tr->current_L = 1;
      }

      // Body-body
      if (i != 0 && fcond < 0.0){
        ri_tr->current_Ks[pindex(i,j,N)] = 1;
        if (ri_tr->encounter_map[i] == 0){
            ri_tr->encounter_map[i] = i;
            ri_tr->encounterN++;
        }
        if (ri_tr->encounter_map[j] == 0){
            ri_tr->encounter_map[j] = j;
            ri_tr->encounterN++;
        }
        if (c0 == 0){
          c0 = 1;
        }
          //if (j<Nactive){ // Two massive particles have a close encounter
          //    ri_tr->tponly_encounter = 0;
          //}
      }
    }
  }

  return c0;
}

int Ftry(struct reb_simulation* const r){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const int Nactive = r->N_active==-1?r->N:r->N_active;
  struct reb_vec3d (*_switch) (const struct reb_simulation* const r, int i, int j) = r->ri_tr.S;
  int ctry = 0;

  for (int i = 0; i < N; i++){
    ri_tr->encounter_map[i] = ri_tr->encounter_map_backup[i];
  }

  for (int i = 0; i < Nactive; i++){
    for (int j = i + 1; j < N; j++){

      struct reb_vec3d cond = _switch(r, i, j);
      double fcond = cond.x;
      double fcond_peri = cond.y;

      // We only care if K_ij goes from 0->1
      // if 1->0 just integrate with BS, probably not worth the hassle
      if (i != 0 && ri_tr->current_Ks[pindex(i,j,N)] == 0 && fcond < 0.0){
        ri_tr->current_Ks[pindex(i,j,N)] = 1;
        if (ri_tr->encounter_map[i]==0){
            ri_tr->encounter_map[i] = i;
            ri_tr->encounterN++;
        }
        if (ri_tr->encounter_map[j]==0){
            ri_tr->encounter_map[j] = j;
            ri_tr->encounterN++;
        }
        ctry = 1;
      }

      // If new pericenter CE has been detected
      if (i == 0 && j < Nactive && ri_tr->current_L == 0 && fcond_peri < 0.0){
        ri_tr->current_L = 1;
        ctry = 1;
      }
    }
  }
  return ctry;
}

int Ftest(struct reb_simulation* const r){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const int Nactive = r->N_active==-1?r->N:r->N_active;
  struct reb_vec3d (*_switch) (const struct reb_simulation* const r, int i, int j) = r->ri_tr.S;
  int ctest = 0;

  for (int i = 1; i < Nactive; i++){
    for (int j = i + 1; j < N; j++){
      struct reb_vec3d cond = _switch(r, i, j);
      double fcond = cond.x;

      // We only care if K_ij goes from 0->1
      // if 1->0 just integrate with BS, probably not worth the hassle
      if (ri_tr->current_Ks[pindex(i,j,N)] == 0 && fcond < 0.0){
        ri_tr->current_Ks[pindex(i,j,N)] = 1;
        if (ctest == 0){
          ctest = 1;
        }
      }
    }
  }
  return ctest;
}

void reb_integrator_trace_select_timestep(struct reb_simulation* const r){
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    double old_dt = r->dt;
    for (int i = ri_tr->nshells - 1; i >= 0; i--){
      if (ri_tr->dt_proposed >= ri_tr->dt_shells[i]){
        r->dt = ri_tr->dt_shells[i];
        ri_tr->current_shell = i;
        break;
      }

      // If we get to the last shell this is the minimum ts
      if (i == 0){
        r->dt = ri_tr->dt_shells[i];
      }
    }
}

// This is Listing 2
void reb_integrator_trace_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    const int N = r->N;
    // Make copy of particles
    memcpy(ri_tr->particles_backup,r->particles,N*sizeof(struct reb_particle));

    // Set encounter map
    ri_tr->encounter_map[0] = 1;
    for (int i=1; i<N; i++){
        ri_tr->encounter_map[i] = 0;
    }
    ri_tr->encounterN = 1;
    ri_tr->current_L = 0;

    for (int i = 0; i < (N-1)*(N-2)/2; i++){
      ri_tr->current_Ks[i] = 0;
      //ri_tr->delta_Ks[i] = 0;
    }

    if (r->t == 0){
      for (int i = 0; i < 6; i++){
        ri_tr->delta_Ks[i] = 0;
      }
    }

    ri_tr->delta_Ks[0] += 1;

    int c0 = F0(r); // Check initial condition
    int ts_rej = 0;

    // TLU TS
    if (c0){
      reb_integrator_trace_select_timestep(r);
    }
    else{
      // no close encounters, we use the lowest timestep for WHFast only
      ri_tr->print = 1;
      r->dt = ri_tr->dt_shells[0];
      ri_tr->dt_proposed = r->dt; // guess for next time we need BS
      //printf("No CE, using lowest timestep %f\n", r->t);
    }



    if (!ri_tr->current_L){ // No pericenter close encounter, proceed as normal

      // Make copy of encounter map if needed
      memcpy(ri_tr->encounter_map_backup,ri_tr->encounter_map,N*sizeof(int));

      // No pericenter CE, so we do the jump step
      if (ri_tr->is_synchronized){
          //printf("synced??\n");
          reb_integrator_trace_jump_step(r, r->dt/2.); // Pdot for B
      }else{
          //printf("not synced\n");
          reb_integrator_trace_jump_step(r, r->dt);
      }
      reb_integrator_trace_interaction_step(r, r->dt/2.);
      reb_integrator_trace_com_step(r,r->dt);

      memcpy(ri_tr->particles_backup_try,r->particles,r->N*sizeof(struct reb_particle));
      reb_integrator_trace_kepler_step(r, r->dt); // We can always advance ALL particles
      ts_rej = reb_trace_bs_step(r, r->dt); // This will do nothing if no close encounters

      // Recheck condition here. Essentially, checking if any Ks changed 0 -> 1
      // continuously recheck condition
      int ctry = Ftry(r);
      memcpy(ri_tr->encounter_map_backup,ri_tr->encounter_map,N*sizeof(int));
      while (ctry || ts_rej){
        //if (ri_tr->print){
        //  printf("Rejected! %f %d %d\n", r->t, ctry, ts_rej);
        //}
        // Either flagged for close encounter or timestep rejection
        // Either way, reset particles
        for (int i=0; i<N; i++){
            r->particles[i] = ri_tr->particles_backup[i];
        }

        if (ts_rej){
          r->dt = ri_tr->dt_shells[ri_tr->current_shell];
        }

        if (ri_tr->is_synchronized){
            reb_integrator_trace_jump_step(r, r->dt/2.); // Pdot for B
        }else{
            reb_integrator_trace_jump_step(r, r->dt);
        }
        reb_integrator_trace_interaction_step(r, r->dt/2.);
        reb_integrator_trace_com_step(r,r->dt);

        memcpy(ri_tr->particles_backup_try,r->particles,r->N*sizeof(struct reb_particle));

        // No Peri CE, WHFast into BS. Add Peri shortly
        reb_integrator_trace_kepler_step(r, r->dt); // We can always advance ALL particles
        ts_rej = reb_trace_bs_step(r, r->dt);
        if (!ts_rej){
          ctry = Ftry(r);
          memcpy(ri_tr->encounter_map_backup,ri_tr->encounter_map,N*sizeof(int));
        }
      }

    }

    else{ // there has been a pericenter close encounter.
      // Immediately integrate entire sim with BS - no Jump step
        //printf("Interaction one\n");
        // BUG: IF YOU START AT PERI CE ERROR IS BAD
      reb_integrator_trace_interaction_step(r, r->dt/2.);
      reb_integrator_trace_com_step(r,r->dt);

      memcpy(ri_tr->particles_backup_try,r->particles,r->N*sizeof(struct reb_particle));

      ri_tr->encounter_map[0] = 1;
      ri_tr->encounterN = N;
      for (int i = 1; i < N; i++){
        ri_tr->encounter_map[i] = i; // Identity map
      }

      ri_tr->encounterNactive = ((r->N_active==-1)?r->N:r->N_active);
      ts_rej = reb_trace_bs_step(r, r->dt); // This will do nothing if no close encounters

      int ctest = Ftest(r);
      if (ctest){
        // Reset to backup values
        for (int i=0; i<N; i++){
            r->particles[i] = ri_tr->particles_backup[i];
        }

        if (ri_tr->is_synchronized){
            reb_integrator_trace_jump_step(r, r->dt/2.); // Pdot for B
        }else{
            reb_integrator_trace_jump_step(r, r->dt);
        }
        reb_integrator_trace_interaction_step(r, r->dt/2.);
        reb_integrator_trace_com_step(r,r->dt);

        memcpy(ri_tr->particles_backup_try,r->particles,r->N*sizeof(struct reb_particle));

        ri_tr->encounterN = N;
        ri_tr->encounter_map[0] = 1; // Identity map except central body
        for (int i = 1; i < N; i++){
          ri_tr->encounter_map[i] = i;
        }
        ri_tr->encounterNactive = ((r->N_active==-1)?r->N:r->N_active);
        ts_rej = reb_trace_bs_step(r, r->dt);
      }
      else{
      }
    }

    reb_integrator_trace_interaction_step(r,r->dt/2.);

    ri_tr->is_synchronized = 0;
    if (ri_tr->safe_mode){
        reb_integrator_trace_synchronize(r);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;

    //printf("Very end of TS: %f %f\n", r->dt, r->t);
}

void reb_integrator_trace_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    if (ri_tr->is_synchronized == 0){
        //printf("synchronizing\n");
        r->gravity = REB_GRAVITY_TRACE; // needed here again for SimulationArchive

        ri_tr->mode=0;
        reb_integrator_trace_jump_step(r,r->dt/2.);
        reb_integrator_trace_dh_to_inertial(r);

        ri_tr->recalculate_coordinates_this_timestep = 1;
        ri_tr->is_synchronized = 1;
    }
}

void reb_integrator_trace_reset(struct reb_simulation* r){
    r->ri_tr.mode = 0;
    r->ri_tr.encounterN = 0;
    r->ri_tr.encounterNactive = 0;
    r->ri_tr.hillfac = 4; // TLu changed to Hernandez (2023)
    r->ri_tr.vSwitch = 0; // TLu changed to Hernandez (2023)
    r->ri_tr.vfac = 3.;
    //r->ri_tr.peri = 0.; // TLu changed to Hernandez (2023)
    r->ri_tr.tponly_encounter = 0;
    r->ri_tr.recalculate_coordinates_this_timestep = 0;
    // Internal arrays (only used within one timestep)
    free(r->ri_tr.particles_backup);
    r->ri_tr.particles_backup = NULL;
    free(r->ri_tr.particles_backup_additionalforces);
    r->ri_tr.particles_backup_additionalforces = NULL;
    free(r->ri_tr.encounter_map);
    r->ri_tr.encounter_map = NULL;
    r->ri_tr.allocatedN = 0;
    r->ri_tr.allocatedN_additionalforces = 0;
    // dcrit array
    free(r->ri_tr.dcrit);
    r->ri_tr.dcrit = NULL;
    r->ri_tr.dcrit_allocatedN = 0;

    free(r->ri_tr.particles_backup_try);
    r->ri_tr.particles_backup_try = NULL;

    free(r->ri_tr.current_Ks);
    r->ri_tr.current_Ks = NULL;
    free(r->ri_tr.delta_Ks);
    r->ri_tr.delta_Ks = NULL;

    r->ri_tr.current_L = 0;

    free(r->ri_tr.encounter_map_backup);
    r->ri_tr.encounter_map_backup = NULL;

    free(r->ri_tr.encounter_map_WH);
    r->ri_tr.encounter_map_WH = NULL;

    free(r->ri_tr.dt_shells);
    r->ri_tr.dt_shells = NULL;

    //free(r->ri_tr.f0);
    //r->ri_tr.f0 = NULL;

    //free(r->ri_tr.f0_peris);
    //r->ri_tr.f0_peris = NULL;

    // free(r->ri_tr.close_encounters);
    // r->close_encounters = NULL;
}

// Check F for pre-timestep using backup particles
/*
int Fp(struct reb_simulation* const r, int* L_rej){
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const double* const dcrit = ri_tr->dcrit;
  const double peri = ri_tr->peri;
  int c0 = 0;

  for (int i = 0; i < N; i++){
    for (int j = i + 1; j < N; j++){
      // pre-timestep
      const double dx0 = ri_tr->particles_backup[i].x - ri_tr->particles_backup[j].x;
      const double dy0 = ri_tr->particles_backup[i].y - ri_tr->particles_backup[j].y;
      const double dz0 = ri_tr->particles_backup[i].z - ri_tr->particles_backup[j].z;
      const double d0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

      double dcritmax = MAX(dcrit[i],dcrit[j]);
      dcritmax *= 1.21;

      double f0 = d0 - dcritmax;
      double f0_peri = d0 - peri;

      // Body-body
      if (i != 0){
        // Populate f0 array so don't need to recalculate
        ri_tr->f0[pindex(i,j,N)] = f0; // probably needs to be fact checked, but seems to be working for now
        if (f0 < 0.0){
          ri_tr->current_Ks[pindex(i,j,N)] = 1;
          c0 = 1;
        }
      }

      // Check for close encounter with central body
      if (i == 0){
          ri_tr->f0_peris[j-1] = f0_peri;
          if (f0_peri < 0.0){
            ri_tr->current_Ls[j-1] = 1;
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
  struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
  const int N = r->N;
  const double* const dcrit = ri_tr->dcrit;
  const double peri = ri_tr->peri;
  int ctry = 0;

  // Ftry only assesses planet-planet
  for (int i = 1; i < N; i++){
    for (int j = i + 1; j < N; j++){
      const double dx = r->particles[i].x - r->particles[j].x;
      const double dy = r->particles[i].y - r->particles[j].y;
      const double dz = r->particles[i].z - r->particles[j].z;
      const double d = sqrt(dx*dx + dy*dy + dz*dz);

      double f0 = ri_tr->f0[pindex(i,j,N)];;
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
*/

// TLu part2 from Hernandez & Dehnen 2023 Listing 3
// simple version: not particle-wise. Integrate whole sim w/ WHFast vs BS
/*
void reb_integrator_trace_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
    const int N = r->N;
    const double* const dcrit = ri_tr->dcrit;

    for (int i = 1; i < N; i++){
      ri_tr->current_Ls[i-1] = 0; // 0 means no close encounter
    }
    ri_tr->current_K = 0;

    // Make copy of particles
    memcpy(ri_tr->particles_backup,r->particles,N*sizeof(struct reb_particle));

    // ytry for each particle
    int L_rej = 0;
    int c0 = Fp(r, &L_rej);

    // If there has been a close encounter, advance whole sim with BS
    if (c0){
      reb_integrator_trace_bs_step(r, r->dt);
    }
    // Otherwise, advance with WHFast
    else{
      reb_integrator_trace_whfast_step(r, r->dt);
    }

    // If close encounter with sun was detected, immediately accept BS step
    // otherwise:
    if (L_rej == 0){ // may be a better way to do this
      // Backup particles from first try in case alternative integrator is rejected
      memcpy(ri_tr->particles_backup_try,r->particles,N*sizeof(struct reb_particle));

      // Now we assess if the step should be accepted by evaluating Ftry
      int ctry = F_cond(r,0);
      if (ctry){
        // Reject first integrator, try alternative integrators
        // In both cases, reset to backup values
        for (int i=0; i<N; i++){
            r->particles[i] = ri_tr->particles_backup[i];
        }

        if (c0){
          // Means we tried BS and rejected. Now try with WHFast
          reb_integrator_trace_whfast_step(r, r->dt);
        }
        else{
          ri_tr->current_K = 1;
          reb_integrator_trace_bs_step(r, r->dt);
        }
        // Now assess Falt
        int calt = F_cond(r,1);
          if (calt){
            // Reject alternative integrator, use original integrators
            for (int i=0; i<N; i++){
                r->particles[i] = ri_tr->particles_backup_try[i];
            }
          }
        }
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}
*/

/*

static void reb_trace_encounter_predict(struct reb_simulation* const r){
  // For now this flags ALL particles TLu
    struct reb_simulation_integrator_trace* ri_tr = &(r->ri_tr);
    struct reb_particle* const particles = r->particles;
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

double reb_integrator_trace_calculate_dcrit_for_particle(struct reb_simulation* r, unsigned int i){
    // Test particles have no dcrit
    if (r->particles[i].m == 0){
      return 0.0;
    }
    struct reb_simulation_integrator_trace* const ri_tr = &(r->ri_tr);
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
    // const double vc = sqrt(GM/fabs(a));
    double dcrit = ri_tr->hillfac*a*cbrt(r->particles[i].m/(3.*r->particles[0].m));
    //printf("%d %f\n", i, dcrit);
    // Criteria 1: average velocity
    //dcrit = MAX(dcrit, vc*0.4*r->dt);
    // Criteria 2: current velocity
    //dcrit = MAX(dcrit, sqrt(v2)*0.4*r->dt);
    // Criteria 3: Hill radius
    // dcrit = MAX(dcrit, ri_tr->hillfac*a*cbrt(r->particles[i].m/(3.*r->particles[0].m)));
    // Criteria 4: physical radius
    //dcrit = MAX(dcrit, 2.*r->particles[i].r);
    return dcrit;
}




*/