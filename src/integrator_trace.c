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
        dcriti = ri_trace->r_crit_hill*di*cbrt(r->particles[i].m/(3.*m0));
    }

    if (r->particles[j].m != 0){

        const double dxj  = r->particles[j].x;  // in dh
        const double dyj  = r->particles[j].y;
        const double dzj  = r->particles[j].z;
        const double dj = sqrt(dxj*dxj + dyj*dyj + dzj*dzj);
        dcritj = ri_trace->r_crit_hill*dj*cbrt(r->particles[j].m/(3.*m0));
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
    const double peri = ri_trace->peri_crit_distance;
    if (peri == 0.0){
      reb_simulation_warning(r,"Pericenter distance parameter not set. Set with r->ri_trace.peri_crit_distance=");
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
    const double pfdot = ri_trace->peri_crit_fdot;
    const double pdist = ri_trace->peri_crit_distance;

    const double dx = r->particles[j].x;
    const double dy = r->particles[j].y;
    const double dz = r->particles[j].z;
    const double d2 = dx*dx + dy*dy + dz*dz;
    const double d = sqrt(d2);

    const double dvx = r->particles[j].vx - r->particles[0].vx;
    const double dvy = r->particles[j].vy - r->particles[0].vy;
    const double dvz = r->particles[j].vz - r->particles[0].vz;

    const double hx = (dy*dvz - dz*dvy);  // specific angular momentum vector
    const double hy = (dz*dvx - dx*dvz);
    const double hz = (dx*dvy - dy*dvx);
    const double h = sqrt ( hx*hx + hy*hy + hz*hz );

    // This only works for bound orbits!
    const double fdot = h / (d2);
    const double peff = (2 * M_PI / fdot); // effective period
    double fcond_peri = peff - pfdot * r->dt;

    // Failsafe: use pericenter pericenter distance
    double fcond_dist = d - pdist;
    return MIN(fcond_peri, fcond_dist);
}

double reb_integrator_trace_peri_switch_none(struct reb_simulation* const r, const unsigned int j){
    // No pericenter flags
    return 1.0;
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
    r->ri_trace.com_pos = com_pos;
    r->ri_trace.com_vel = com_vel;
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
    const int current_C = ri_trace->current_C;
    if (current_C) return; // No jump step for pericenter approaches

    const int N_active = r->N_active==-1?r->N:r->N_active;

    // If TP type 1, use r->N. Else, use N_active.
    const int N = r->testparticle_type==0 ? N_active: r->N;

    double px=0., py=0., pz=0.;
    for (int i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m;
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px *= dt/r->particles[0].m;
    py *= dt/r->particles[0].m;
    pz *= dt/r->particles[0].m;

    const int N_all = r->N;
    for (int i=1;i<N_all;i++){
        particles[i].x += px;
        particles[i].y += py;
        particles[i].z += pz;
    }
}

void reb_integrator_trace_com_step(struct reb_simulation* const r, double dt){
    r->ri_trace.com_pos.x += dt*r->ri_trace.com_vel.x;
    r->ri_trace.com_pos.y += dt*r->ri_trace.com_vel.y;
    r->ri_trace.com_pos.z += dt*r->ri_trace.com_vel.z;
}

void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    //struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        reb_whfast_kepler_solver(r,r->particles,r->G*r->particles[0].m,i,dt);
    }
}

void reb_integrator_trace_update_particles(struct reb_simulation* r, const double* y){
    int N = r->ri_trace.encounter_N;
    int* map = r->ri_trace.encounter_map;

    for (int i=0; i<N; i++){
        int mi = map[i];
        struct reb_particle* const p = &(r->particles[mi]);
        p->x  = y[i*6+0];
        p->y  = y[i*6+1];
        p->z  = y[i*6+2];
        p->vx = y[i*6+3];
        p->vy = y[i*6+4];
        p->vz = y[i*6+5];
    }
}

void reb_integrator_trace_nbody_derivatives(struct reb_ode* ode, double* const yDot, const double* const y, double const t){
    struct reb_simulation* const r = ode->r;
    // TRACE always needs this to ensure the right Hamiltonian is evolved
    reb_integrator_trace_update_particles(r, y);
    reb_simulation_update_acceleration(r);

    // TLu Levison & Duncan 22, 23 EoMs
    double px=0., py=0., pz=0.;
    int* map = r->ri_trace.encounter_map;
    int N = r->ri_trace.encounter_N;

    if (map==NULL){
        reb_simulation_error(r, "Cannot access TRACE map from BS.");
        return;
    }

    // Kepler Step
    // This is only for pericenter approach
    if (r->ri_trace.current_C){
        for (int i=1;i<r->N;i++){ // all particles
            px += r->particles[i].vx*r->particles[i].m; // in dh
            py += r->particles[i].vy*r->particles[i].m;
            pz += r->particles[i].vz*r->particles[i].m;
        }
        px /= r->particles[0].m;
        py /= r->particles[0].m;
        pz /= r->particles[0].m;

    }
    yDot[0*6+0] = 0.0;
    yDot[0*6+1] = 0.0;
    yDot[0*6+2] = 0.0;
    yDot[0*6+3] = 0.0;
    yDot[0*6+4] = 0.0;
    yDot[0*6+5] = 0.0;

    for (int i=1; i<N; i++){
        int mi = map[i];
        const struct reb_particle p = r->particles[mi];
        yDot[i*6+0] = p.vx + px; // Already checked for current_L
        yDot[i*6+1] = p.vy + py;
        yDot[i*6+2] = p.vz + pz;
        yDot[i*6+3] = p.ax;
        yDot[i*6+4] = p.ay;
        yDot[i*6+5] = p.az;
    }
}

void reb_integrator_trace_bs_step(struct reb_simulation* const r, double dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);

    if (ri_trace->encounter_N < 2){
        // No close encounters, skip
        return;
    }

    int i_enc = 0;
    ri_trace->encounter_N_active = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(ri_trace->encounter_map_internal[i]){
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = ri_trace->particles_backup_kepler[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            ri_trace->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
                ri_trace->encounter_N_active++;
                if (ri_trace->tponly_encounter){
                    ri_trace->particles_backup_kepler[i] = tmp;         // Make copy of particles after the kepler step.
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
    const double t_needed = r->t + dt;
    reb_integrator_bs_reset(r);

    struct reb_ode* nbody_ode = reb_ode_create(r, ri_trace->encounter_N*3*2);
    nbody_ode->derivatives = reb_integrator_trace_nbody_derivatives;
    nbody_ode->needs_nbody = 0;

    // TODO: Support backwards integrations
    while(r->t < t_needed && fabs(dt/old_dt)>1e-14 ){
        double* y = nbody_ode->y;
        
        // In case of overshoot
        if (r->t + dt >  t_needed){
            dt = t_needed - r->t;
        }

        struct reb_particle star = r->particles[0]; // backup velocity
        r->particles[0].vx = 0; // star does not move in dh
        r->particles[0].vy = 0;
        r->particles[0].vz = 0;

        for (unsigned int i=0; i<ri_trace->encounter_N; i++){
            const int mi = ri_trace->encounter_map[i];
            const struct reb_particle p = r->particles[mi];
            y[i*6+0] = p.x;
            y[i*6+1] = p.y;
            y[i*6+2] = p.z;
            y[i*6+3] = p.vx;
            y[i*6+4] = p.vy;
            y[i*6+5] = p.vz;
        }

        int success = reb_integrator_bs_step(r, dt);
        if (success){
            r->t += dt;
        }
        dt = r->ri_bs.dt_proposed;
        reb_integrator_trace_update_particles(r, nbody_ode->y);

        r->particles[0].vx = star.vx; // restore every timestep for collisions
        r->particles[0].vy = star.vy;
        r->particles[0].vz = star.vz;

        reb_collision_search(r);

        if (nbody_ode->length != ri_trace->encounter_N*3*2){
            if (ri_trace->encounter_N*3*2 > nbody_ode->N_allocated){
                reb_simulation_error(r, "Cannot add particles during encounter step");
            }
            nbody_ode->length = ri_trace->encounter_N*3*2;
            r->ri_bs.first_or_last_step = 1;
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

    // if only test particles encountered massive bodies, reset the
    // massive body coordinates to their post Kepler step state
    if(ri_trace->tponly_encounter){
        for (unsigned int i=1; i < ri_trace->encounter_N_active; i++){
            unsigned int mi = ri_trace->encounter_map[i];
            r->particles[mi] = ri_trace->particles_backup_kepler[mi];
        }
    }

    reb_ode_free(nbody_ode);

    r->t = old_t;
    ri_trace->mode = 0;
}

void reb_integrator_trace_kepler_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    memcpy(ri_trace->particles_backup_kepler,r->particles,r->N*sizeof(struct reb_particle));
    reb_integrator_trace_whfast_step(r, _dt);
    reb_integrator_trace_bs_step(r, _dt);
}


void reb_integrator_trace_part1(struct reb_simulation* r){
    if (r->N_var_config){
        reb_simulation_warning(r,"TRACE does not work with variational equations.");
    }

    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;

    if (ri_trace->N_allocated<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_trace->particles_backup       = realloc(ri_trace->particles_backup,sizeof(struct reb_particle)*N);

        if (ri_trace->current_Ks){ // Free previously allocated matrix
            for (int k = 0; k < ri_trace->N_allocated; ++k) {
                free(ri_trace->current_Ks[k]);
            }
            free(ri_trace->current_Ks);
        }
        ri_trace->current_Ks             = malloc(sizeof(int*)*N); // This is inefficient for now, can be Nactive instead of N
        for (int k = 0; k < N; ++k) {
            ri_trace->current_Ks[k]      = malloc(sizeof(int)*N);
        }

        ri_trace->encounter_map          = realloc(ri_trace->encounter_map,sizeof(int)*N);
        ri_trace->encounter_map_internal = realloc(ri_trace->encounter_map_internal,sizeof(int)*N); // Do we need this now?

        // Only need this stuff for Listing 3
        ri_trace->particles_backup_kepler   = realloc(ri_trace->particles_backup_kepler,sizeof(struct reb_particle)*N);
        ri_trace->N_allocated = N;
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

    r->gravity = REB_GRAVITY_TRACE;
    ri_trace->mode = 2; // Do not calculate gravity in-between timesteps. TRACE will call reb_update_acceleration itself.

    reb_integrator_trace_inertial_to_dh(r);

}

void reb_integrator_trace_pre_ts_check(struct reb_simulation* const r){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;
    double (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = r->ri_trace.S;
    double (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = r->ri_trace.S_peri;
    
    // Clear encounter maps
    for (unsigned int i=0; i<r->N; i++){
        ri_trace->encounter_map[i] = 0;
        ri_trace->encounter_map_internal[i] = 0;
    }
    ri_trace->encounter_map_internal[0] = 1;

    // Reset encounter triggers.
    ri_trace->encounter_N = 1;
    ri_trace->current_C = 0;

    for (int i = 0; i < N; i++){
        for (unsigned int j = i + 1; j < N; j++){
            ri_trace->current_Ks[i][j] = 0;
        }
    }

    if (r->testparticle_type == 1){
        ri_trace->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        ri_trace->tponly_encounter = 1;
    }

    // Check for pericenter CE
    for (int j = 1; j < Nactive; j++){
        if (_switch_peri(r, j) < 0.0){
            ri_trace->current_C = 1;
            if (j < Nactive){ // Two massive particles have a close encounter
                ri_trace->tponly_encounter = 0;
                break; // No need to check other particles
            }
        }
    }
    
    if (ri_trace->current_C){
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        ri_trace->encounter_N = N;
        for (int i = 0; i < N; i++){
            ri_trace->encounter_map_internal[i] = 1; //  trigger encounter
        }
        ri_trace->encounter_N_active = ((r->N_active==-1)?r->N:r->N_active);
        return; // No need to check other condition
    }

    // Body-body
    // there cannot be TP-TP CEs
    for (int i = 1; i < Nactive; i++){
        for (int j = i + 1; j < N; j++){
            if (_switch(r, i, j) < 0.0){
                ri_trace->current_Ks[i][j] = 1;
                if (ri_trace->encounter_map_internal[i] == 0){
                    ri_trace->encounter_map_internal[i] = 1; // trigger encounter
                    ri_trace->encounter_N++;
                }
                if (ri_trace->encounter_map_internal[j] == 0){
                    ri_trace->encounter_map_internal[j] = 1; // trigger encounter
                    ri_trace->encounter_N++;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    ri_trace->tponly_encounter = 0;
                }
            }
        }
    }
}

double reb_integrator_trace_post_ts_check(struct reb_simulation* const r){
    // This function returns 1 if any new encounters occured.
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;
    double (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = r->ri_trace.S;
    double (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = r->ri_trace.S_peri;
    int new_close_encounter = 0; // New CEs

    // Check for pericenter CE
    for (int j = 1; j < Nactive; j++){
        if (_switch_peri(r, j) < 0.0){
            ri_trace->current_C = 1;
            new_close_encounter = 1;

            if (j < Nactive){ // Two massive particles have a close encounter
                ri_trace->tponly_encounter = 0;
                break; // No need to check other particles
            }
        }
    }
    if (ri_trace->current_C){
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        ri_trace->encounter_N = N;
        for (int i = 0; i < N; i++){
            ri_trace->encounter_map_internal[i] = 1; // trigger encounter
        }
        ri_trace->encounter_N_active = ((r->N_active==-1)?r->N:r->N_active);
        return new_close_encounter; // No need to check for particle particle encounters
    }


    // Body-body
    // there cannot be TP-TP CEs
    for (int i = 1; i < Nactive; i++){
        for (int j = i + 1; j < N; j++){
            if (ri_trace->current_Ks[i][j] == 0){
              if (_switch(r, i, j) < 0.0){
                  ri_trace->current_Ks[i][j] = 1;
                  new_close_encounter = 1;
                  if (ri_trace->encounter_map_internal[i] == 0){
                      ri_trace->encounter_map_internal[i] = 1; // trigger encounter
                      ri_trace->encounter_N++;
                  }
                  if (ri_trace->encounter_map_internal[j] == 0){
                      ri_trace->encounter_map_internal[j] = 1; // trigger encounter
                      ri_trace->encounter_N++;
                  }

                  if (j < Nactive){ // Two massive particles have a close encounter
                      ri_trace->tponly_encounter = 0;
                  }
              }
            }
        }
    }
    return new_close_encounter;
}

void reb_integrator_trace_part2(struct reb_simulation* const r){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    
    ri_trace->mode = 0;
                        
    // This will be set to 1 if a collision occured.
    ri_trace->force_accept = 0;

    // Create copy of all particle to allow for the step to be rejected.
    memcpy(ri_trace->particles_backup, r->particles, N*sizeof(struct reb_particle));

    // Check if there are any close encounters
    reb_integrator_trace_pre_ts_check(r);

    reb_integrator_trace_interaction_step(r, r->dt/2.);
    reb_integrator_trace_jump_step(r, r->dt/2.);
    reb_integrator_trace_kepler_step(r, r->dt);
    reb_integrator_trace_com_step(r,r->dt);
    reb_integrator_trace_jump_step(r, r->dt/2.);
    reb_integrator_trace_interaction_step(r, r->dt/2.);

    // We might need to check again for close encounters to ensure time reversibility. However:
    // - We alaways accept the step if a collision occured
    //   as it is impossible to undo the collision.
    // - We don't need to check the encounter conditions if the pericenter 
    //   condition was already triggered because all particles are integrated with BS.
    if (!ri_trace->force_accept && !ri_trace->current_C){
        // Check for new close encounters not present pre timestep
        if (reb_integrator_trace_post_ts_check(r)){
            // New encounters were found. Will reject the step.
            // Revert particles to the beginning of the step.
            memcpy(r->particles, ri_trace->particles_backup, N*sizeof(struct reb_particle));

            // Do step again
            reb_integrator_trace_interaction_step(r, r->dt/2.);
            reb_integrator_trace_jump_step(r, r->dt/2.);
            reb_integrator_trace_kepler_step(r, r->dt);
            reb_integrator_trace_com_step(r,r->dt);
            reb_integrator_trace_jump_step(r, r->dt/2.);
            reb_integrator_trace_interaction_step(r, r->dt/2.);
        }
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;

    reb_integrator_trace_dh_to_inertial(r);
}

void reb_integrator_trace_synchronize(struct reb_simulation* r){
}

void reb_integrator_trace_reset(struct reb_simulation* r){
    r->ri_trace.mode = 0;
    r->ri_trace.encounter_N = 0;
    r->ri_trace.encounter_N_active = 0;
    r->ri_trace.r_crit_hill = 3;
    r->ri_trace.peri_crit_fdot = 17.;
    r->ri_trace.peri_crit_distance = 0.;

    // Internal arrays (only used within one timestep)
    free(r->ri_trace.particles_backup);
    r->ri_trace.particles_backup = NULL;
    free(r->ri_trace.particles_backup_kepler);
    r->ri_trace.particles_backup_kepler = NULL;
    free(r->ri_trace.particles_backup_additional_forces);
    r->ri_trace.particles_backup_additional_forces = NULL;


    free(r->ri_trace.encounter_map);
    r->ri_trace.encounter_map = NULL;
    free(r->ri_trace.encounter_map_internal);
    r->ri_trace.encounter_map_internal = NULL;

    r->ri_trace.current_C = 0;
    if (r->ri_trace.current_Ks){
        for (int k=0; k < r->ri_trace.N_allocated; k++) {
            free(r->ri_trace.current_Ks[k]);
        }
        free(r->ri_trace.current_Ks);
        r->ri_trace.current_Ks = NULL;
    }

    r->ri_trace.S = NULL;
    r->ri_trace.S_peri = NULL;
    
    r->ri_trace.N_allocated = 0;
    r->ri_trace.N_allocated_additional_forces = 0;
}
