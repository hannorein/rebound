/**
 * @file    integrator_trace.c
 * @brief   TRACE
 * @author  Tiger Lu
 *
 * @section LICENSE
 * Copyright (c) 2023 Tiger Lu
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

double reb_integrator_trace_switch_default(struct reb_simulation* const r, const unsigned int i, const unsigned int j){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    double dcriti = 0.0;
    double dcritj = 0.0;

    const double m0 = r->particles[0].m;

    if (r->particles[i].m != 0){

        const double dxi  = r->particles[i].x;  // in dh
        const double dyi  = r->particles[i].y;
        const double dzi  = r->particles[i].z;
        const double di = sqrt(dxi*dxi + dyi*dyi + dzi*dzi);
        dcriti = ri_trace->hillfac*di*cbrt(r->particles[i].m/(3.*m0));
    }

    if (r->particles[j].m != 0){

        const double dxj  = r->particles[j].x;  // in dh
        const double dyj  = r->particles[j].y;
        const double dzj  = r->particles[j].z;
        const double dj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);
        dcritj = ri_trace->hillfac*dj*cbrt(r->particles[j].m/(3.*m0));
    }

    const double dx = r->particles[i].x - r->particles[j].x;
    const double dy = r->particles[i].y - r->particles[j].y;
    const double dz = r->particles[i].z - r->particles[j].z;
    const double d = sqrt(dx*dx + dy*dy + dz*dz);


    // Use traditional switching function
    double dcritmax = MAX(dcriti,dcritj);
    dcritmax *= 1.21;

    double fcond = d - dcritmax;
    return fcond;
}

double reb_integrator_trace_peri_switch_distance(struct reb_simulation* const r, const unsigned int j){
    const struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const double peri = ri_trace->peri_distance;
    if (peri == 0.0){
      reb_simulation_warning(r,"Pericenter condition not set. Set with r->ri_trace.peri_distance=");
    }

    const double dx = r->particles[0].x - r->particles[j].x;
    const double dy = r->particles[0].y - r->particles[j].y;
    const double dz = r->particles[0].z - r->particles[j].z;
    const double d = sqrt(dx*dx + dy*dy + dz*dz);

    double fcond_peri = d - peri;
    return fcond_peri;
}


double reb_integrator_trace_peri_switch_default(struct reb_simulation* const r, const unsigned int j){
    const struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const double pfdot = ri_trace->peri_fdot;
    const double pdist = ri_trace->peri_distance;

    const double dx = r->particles[j].x;
    const double dy = r->particles[j].y;
    const double dz = r->particles[j].z;
    const double d2 = dx*dx + dy*dy + dz*dz;
    const double d = sqrt(d2);

    const double dvx = r->particles[j].vx - r->particles[0].vx;
    const double dvy = r->particles[j].vy - r->particles[0].vy;
    const double dvz = r->particles[j].vz - r->particles[0].vz;
    const double v2 = dvx * dvx + dvy * dvy + dvz * dvz;

    const double hx = (dy*dvz - dz*dvy);  // specific angular momentum vector
    const double hy = (dz*dvx - dx*dvz);
    const double hz = (dx*dvy - dy*dvx);
    const double h = sqrt ( hx*hx + hy*hy + hz*hz );

    // This only works for bound orbits!
    const double fdot = h / (d2);
    const double peff = (2 * M_PI / fdot); // effective period
    double fcond_peri = (peff / r->dt) - pfdot;

    // Failsafe - use velocity dependent condition
    //double f_vel = d / sqrt(3. * v2 + (r->G * (r->particles[0].m + r->particles[j].m) / d));
    //double fcond_vel = (f_vel/r->dt) - 1.;

    // Failsafe: use pericenter pericenter distance
    double fcond_dist = 0.0;
    if (pdist > 0.0){
        fcond_dist = d - pdist;
    }

    return MIN(fcond_peri, fcond_dist);
}

void reb_integrator_trace_inertial_to_dh(struct reb_simulation* r){
    //printf("Swapped to DH at %f\n", r->t);
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
    r->ri_trace.com_pos = com_pos;
    r->ri_trace.com_vel = com_vel;
}

void reb_integrator_trace_dh_to_inertial(struct reb_simulation* r){
    //printf("Swapped to inertial at %f\n", r->t);
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
    particles[0].x = r->ri_trace.com_pos.x - temp.x;
    particles[0].y = r->ri_trace.com_pos.y - temp.y;
    particles[0].z = r->ri_trace.com_pos.z - temp.z;

    for (int i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += r->ri_trace.com_vel.x;
        particles[i].vy += r->ri_trace.com_vel.y;
        particles[i].vz += r->ri_trace.com_vel.z;
    }
    particles[0].vx = r->ri_trace.com_vel.x - temp.vx;
    particles[0].vy = r->ri_trace.com_vel.y - temp.vy;
    particles[0].vz = r->ri_trace.com_vel.z - temp.vz;
}

void reb_integrator_trace_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    reb_simulation_update_acceleration(r);
    for (int i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_trace_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;

    struct reb_integrator_trace* ri_trace = &(r->ri_trace);
    const int current_L = ri_trace->current_L;

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
    r->ri_trace.com_pos.x += dt*r->ri_trace.com_vel.x;
    r->ri_trace.com_pos.y += dt*r->ri_trace.com_vel.y;
    r->ri_trace.com_pos.z += dt*r->ri_trace.com_vel.z;
}

// Old Kepler
void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    //struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        //if (ri_trace->encounter_map[i] != 0){
        reb_whfast_kepler_solver(r,r->particles,r->G*r->particles[0].m,i,dt); // in dh
                                                                              //}
    }
}

void reb_integrator_trace_bs_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);

    if (ri_trace->encounterN < 2){
        // No close encounters, skip
        return;
    }

    int i_enc = 0;
    ri_trace->encounterNactive = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(ri_trace->encounter_map_internal[i]){
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = ri_trace->particles_backup_try[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            ri_trace->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
                ri_trace->encounterNactive++;
                if (ri_trace->tponly_encounter){
                    ri_trace->particles_backup_try[i] = tmp;         // Make copy of particles after the kepler step.
                                                                  // used to restore the massive objects' states in the case
                                                                  // of only massless test-particle encounters
                }
            }
        }
    }

    ri_trace->mode = 1;
    // run
    const double old_dt = r->dt;
    const double old_t = r->t;
    double t_needed = r->t + _dt;
    //reb_integrator_bs_reset(r);

    r->dt = _dt; // start with a small timestep.

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
    if(ri_trace->tponly_encounter){
        for (unsigned int i=1; i < ri_trace->encounterNactive; i++){
            unsigned int mi = ri_trace->encounter_map[i];
            r->particles[mi] = ri_trace->particles_backup_try[mi];
        }
    }

    r->t = old_t;
    r->dt = old_dt;
    ri_trace->mode = 0;
    // return reject;
}

void reb_integrator_trace_kepler_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    memcpy(ri_trace->particles_backup_try,r->particles,r->N*sizeof(struct reb_particle));
    reb_integrator_trace_whfast_step(r, _dt);
    reb_integrator_trace_bs_step(r, _dt);
    // return rej;
}


void reb_integrator_trace_part1(struct reb_simulation* r){
    //printf("Start: %e %e %e %e %e %e\n", r->particles[1].x,r->particles[1].y,r->particles[1].z,r->particles[1].vx,r->particles[1].vy,r->particles[1].vz);
    if (r->N_var_config){
        reb_simulation_warning(r,"TRACE does not work with variational equations.");
    }

    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;

    if (ri_trace->allocatedN<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_trace->particles_backup       = realloc(ri_trace->particles_backup,sizeof(struct reb_particle)*N);
        ri_trace->current_Ks             = realloc(ri_trace->current_Ks, sizeof(int*)*N); // This is inefficient for now, can be Nactive instead of N
        for (int k = 0; k < N; ++k) {
            ri_trace->current_Ks[k]      = realloc(ri_trace->current_Ks[k], sizeof(int)*N);
        }

        ri_trace->encounter_map          = realloc(ri_trace->encounter_map,sizeof(int)*N);
        ri_trace->encounter_map_internal = realloc(ri_trace->encounter_map_internal,sizeof(int)*N); // Do we need this now?

        // Only need this stuff for Listing 3
        ri_trace->particles_backup_try   = realloc(ri_trace->particles_backup_try,sizeof(struct reb_particle)*N);
        ri_trace->allocatedN = N;
    }

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_simulation_warning(r,"TRACE only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_TRACE){
        reb_simulation_warning(r,"TRACE has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }

    if (ri_trace->S == NULL){
        ri_trace->S = reb_integrator_trace_switch_default;
    }

    if (ri_trace->S_peri == NULL){
        ri_trace->S_peri = reb_integrator_trace_peri_switch_default;
    }

    // Set nbody ODE here to know if we need to integrate other ODEs in step one.
    if (r->ri_bs.nbody_ode == NULL){
        r->ri_bs.nbody_index = r->N_odes;
        r->ri_bs.nbody_ode = reb_ode_create(r, r->N*3*2);
        r->ri_bs.nbody_ode->derivatives = reb_integrator_bs_nbody_derivatives;
        r->ri_bs.nbody_ode->needs_nbody = 0; // No need to update unless there's another ode
        r->ri_bs.first_or_last_step = 1;
    }

    // I can't think of a better way to do this, but TRACE needs to kill needs_nbody for each additional ODE
    if (r->N_odes > 1){
      for (int s=0; s < r->N_odes; s++){
        r->odes[s]->needs_nbody = 0;
      }
    }


    r->gravity = REB_GRAVITY_TRACE;
    ri_trace->mode = 0;
    ri_trace->force_accept = 0;

    //printf("Start: %e %e %e %e %e %e\n", r->particles[1].x,r->particles[1].y,r->particles[1].z,r->particles[1].vx,r->particles[1].vy,r->particles[1].vz);
    reb_integrator_trace_inertial_to_dh(r);
    //printf("Post Transform TO: %e %e %e %e %e %e\n", r->particles[1].x,r->particles[1].y,r->particles[1].z,r->particles[1].vx,r->particles[1].vy,r->particles[1].vz);

    // Clear encounter maps
    for (unsigned int i=0; i<r->N; i++){
        ri_trace->encounter_map[i] = 0;
        ri_trace->encounter_map_internal[i] = 0;
    }
    ri_trace->encounter_map_internal[0] = 1;
}

int reb_integrator_trace_Fcond(struct reb_simulation* const r){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;

    int new_c = 0; // New CEs

    // Switching functions
    double (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = r->ri_trace.S;
    double (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = r->ri_trace.S_peri;

    if (r->testparticle_type == 1){
        ri_trace->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        ri_trace->tponly_encounter = 1;
    }

    // Check for pericenter CE
    // test particles cannot have pericenter CEs OR CAN THEY?? Need to resolve pericenter
    for (int j = 1; j < Nactive; j++){
        double fcond_peri = _switch_peri(r, j);
        if (fcond_peri < 0.0 && ri_trace->current_L == 0){
            ri_trace->current_L = 1;
            new_c = 1;

            if (j < Nactive){ // Two massive particles have a close encounter
                ri_trace->tponly_encounter = 0;
            }
        }
    }
    //exit(1);

    // Body-body
    // there cannot be TP-TP CEs
    for (int i = 1; i < Nactive; i++){
        for (int j = i + 1; j < N; j++){

            double fcond = _switch(r, i, j);

            if (fcond < 0.0){
                if (ri_trace->encounter_map_internal[i] == 0){
                    ri_trace->encounter_map_internal[i] = i;
                    ri_trace->encounterN++;
                }
                if (ri_trace->encounter_map_internal[j] == 0){
                    ri_trace->encounter_map_internal[j] = j;
                    ri_trace->encounterN++;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    ri_trace->tponly_encounter = 0;
                }

                // Checks for switching Kij 0->1. Initialized as all 0 the first time of asking.
                if (ri_trace->current_Ks[i][j] == 0){
                    ri_trace->current_Ks[i][j] = 1;
                    new_c = 1;
                    //if (ri_trace->print){
                    //printf("New CE at %f detected between %d %d\n", r->t,i, j);
                    //exit(1);
                    //}
                }
            }
        }
    }

    return new_c;
}

// This is Listing 2
void reb_integrator_trace_part2(struct reb_simulation* const r){
    //printf("TRACE part 2\n");
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    // Make copy of particles
    memcpy(ri_trace->particles_backup,r->particles,N*sizeof(struct reb_particle));
    ri_trace->encounterN = 1;
    ri_trace->current_L = 0;

    for (int i = 0; i < N; i++){
        for (unsigned int j = i + 1; j < N; j++){
            ri_trace->current_Ks[i][j] = 0;
        }
    }

    reb_integrator_trace_Fcond(r); // return value ignored.

    if (ri_trace->current_L){ //more efficient way to check if we need to redo this...
                           // Pericenter close encounter detected. We integrate the entire simulation with BS
        //printf("\nDetected peri CE\n");
        ri_trace->encounter_map_internal[0] = 1;
        ri_trace->encounterN = N;
        for (int i = 1; i < N; i++){
            ri_trace->encounter_map_internal[i] = i; // Identity map
        }
        ri_trace->encounterNactive = ((r->N_active==-1)?r->N:r->N_active);
    }

    reb_integrator_trace_interaction_step(r, r->dt/2.);
    reb_integrator_trace_jump_step(r, r->dt/2.);
    reb_integrator_trace_kepler_step(r, r->dt); // always accept this
    reb_integrator_trace_com_step(r,r->dt);
    reb_integrator_trace_jump_step(r, r->dt/2.);
    reb_integrator_trace_interaction_step(r, r->dt/2.);

    // Check for new close_encounters
    if (reb_integrator_trace_Fcond(r) && !ri_trace->force_accept){
        // REJECT STEP
        // reset simulation and try again with new timestep
        for (int i=0; i<N; i++){
            // Reject & reset
            r->particles[i] = ri_trace->particles_backup[i];
        }

        if (ri_trace->current_L){ //more efficient way to check if we need to redo this...
                               // Pericenter close encounter detected. We integrate the entire simulation with BS
            ri_trace->encounter_map_internal[0] = 1;
            ri_trace->encounterN = N;
            for (int i = 1; i < N; i++){
                ri_trace->encounter_map_internal[i] = i; // Identity map
            }
            ri_trace->encounterNactive = ((r->N_active==-1)?r->N:r->N_active);
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
    ri_trace->mode = 0;

    r->gravity = REB_GRAVITY_TRACE; // Is this needed?

    reb_integrator_trace_dh_to_inertial(r);
}

void reb_integrator_trace_synchronize(struct reb_simulation* r){
}

void reb_integrator_trace_reset(struct reb_simulation* r){
    r->ri_trace.mode = 0;
    r->ri_trace.encounterN = 0;
    r->ri_trace.encounterNactive = 0;
    r->ri_trace.hillfac = 4; // TLu changed to Hernandez (2023)
    r->ri_trace.peri_fdot = 16.;
    r->ri_trace.peri_distance = 0.;

    //r->ri_trace.peri = 0.; // TLu changed to Hernandez (2023)
    // Internal arrays (only used within one timestep)
    free(r->ri_trace.particles_backup);
    r->ri_trace.particles_backup = NULL;

    free(r->ri_trace.encounter_map);
    r->ri_trace.encounter_map = NULL;
    r->ri_trace.allocatedN = 0;
    r->ri_trace.allocatedN_additionalforces = 0;

    free(r->ri_trace.particles_backup_try);
    r->ri_trace.particles_backup_try = NULL;

    if (r->ri_trace.current_Ks){
        for (int k = 0; k < r->N; ++k) {
            r->ri_trace.current_Ks[k] = NULL;
        }
        free(r->ri_trace.current_Ks);
        r->ri_trace.current_Ks = NULL;
    }

    r->ri_trace.current_L = 0;

    free(r->ri_trace.encounter_map_internal);
    r->ri_trace.encounter_map_internal = NULL;

    r->ri_trace.S = NULL;
    r->ri_trace.S_peri = NULL;

}
