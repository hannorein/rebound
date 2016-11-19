/**
 * @file    integrator_whfasthelio.c
 * @brief   WHFASTHELIO integration scheme.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the WHFast integration scheme in
 *          Heliocentric Coordinates.
 *          Based on WHFast, described in Rein & Tamayo 2015.
 * 
 * @section LICENSE
 * Copyright (c) 2016 Hanno Rein, Daniel Tamayo
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
#include <sys/time.h>
#include "rebound.h"
#include "particle.h"
#include "gravity.h"
#include "integrator.h"
#include "integrator_whfast.h"
#include "integrator_whfasthelio.h"

/****************************** 
 * Coordinate transformations */
static void to_helio_posvel(const struct reb_particle* const particles, struct reb_particle* const p_h, const int N){
    p_h[0].x  = 0.;
    p_h[0].y  = 0.;
    p_h[0].z  = 0.;
    p_h[0].vx = 0.;
    p_h[0].vy = 0.;
    p_h[0].vz = 0.;
    p_h[0].m  = 0.;
    for (unsigned int i=0;i<N;i++){
        p_h[0].x  += particles[i].x *particles[i].m;
        p_h[0].y  += particles[i].y *particles[i].m;
        p_h[0].z  += particles[i].z *particles[i].m;
        p_h[0].vx += particles[i].vx*particles[i].m;
        p_h[0].vy += particles[i].vy*particles[i].m;
        p_h[0].vz += particles[i].vz*particles[i].m;
        p_h[0].m  += particles[i].m;
    }
    p_h[0].x  /= p_h[0].m;
    p_h[0].y  /= p_h[0].m;
    p_h[0].z  /= p_h[0].m;
    p_h[0].vx /= p_h[0].m;
    p_h[0].vy /= p_h[0].m;
    p_h[0].vz /= p_h[0].m;

    for (unsigned int i=1;i<N;i++){
        p_h[i].x  = particles[i].x  - particles[0].x ; 
        p_h[i].y  = particles[i].y  - particles[0].y ;
        p_h[i].z  = particles[i].z  - particles[0].z ;
        p_h[i].vx = particles[i].vx - p_h[0].vx;
        p_h[i].vy = particles[i].vy - p_h[0].vy;
        p_h[i].vz = particles[i].vz - p_h[0].vz;
        p_h[i].m  = particles[i].m;
    }
}

static void to_inertial_pos(struct reb_particle* const particles, const struct reb_particle* const p_h, const int N){
    const double mtot = p_h[0].m;
    particles[0].x  = p_h[0].x;
    particles[0].y  = p_h[0].y;
    particles[0].z  = p_h[0].z;
    for (unsigned int i=1;i<N;i++){
        particles[0].x  -= p_h[i].x*particles[i].m/mtot;
        particles[0].y  -= p_h[i].y*particles[i].m/mtot;
        particles[0].z  -= p_h[i].z*particles[i].m/mtot;
    }
    for (unsigned int i=1;i<N;i++){
        particles[i].x = p_h[i].x+particles[0].x;
        particles[i].y = p_h[i].y+particles[0].y;
        particles[i].z = p_h[i].z+particles[0].z;
    }
}
static void to_inertial_posvel(struct reb_particle* const particles, const struct reb_particle* const p_h, const int N){
    to_inertial_pos(particles,p_h,N);
    const double mtot = p_h[0].m;
    const double m0 = particles[0].m;
    for (unsigned int i=1;i<N;i++){
        particles[i].vx = p_h[i].vx+p_h[0].vx;
        particles[i].vy = p_h[i].vy+p_h[0].vy;
        particles[i].vz = p_h[i].vz+p_h[0].vz;
    }
    particles[0].vx = p_h[0].vx*mtot/m0;
    particles[0].vy = p_h[0].vy*mtot/m0;
    particles[0].vz = p_h[0].vz*mtot/m0;
    for (unsigned int i=1;i<N;i++){
        particles[0].vx -= particles[i].vx*particles[i].m/m0;
        particles[0].vy -= particles[i].vy*particles[i].m/m0;
        particles[0].vz -= particles[i].vz*particles[i].m/m0;
    }
}

/***************************** 
 * Operators                 */
static void reb_whfasthelio_jump_step(const struct reb_simulation* const r, double _dt){
    const double m0 = r->particles[0].m;
    const int N_real = r->N-r->N_var;
    struct reb_particle* const p_h = r->ri_whfasthelio.p_h;
    double px = 0;
    double py = 0;
    double pz = 0;
    for (unsigned int i=1;i<N_real;i++){
        px += p_h[i].m* p_h[i].vx;
        py += p_h[i].m* p_h[i].vy;
        pz += p_h[i].m* p_h[i].vz;
    }
    for (unsigned int i=1;i<N_real;i++){
        p_h[i].x += _dt * px/m0;
        p_h[i].y += _dt * py/m0;
        p_h[i].z += _dt * pz/m0;
    }
}

static void reb_whfasthelio_interaction_step(const struct reb_simulation* const r, const double _dt){
    struct reb_particle* particles = r->particles;
    const int N_real = r->N-r->N_var;
    struct reb_particle* const p_h = r->ri_whfasthelio.p_h;
    for (unsigned int i=1;i<N_real;i++){
        p_h[i].vx += _dt*particles[i].ax;
        p_h[i].vy += _dt*particles[i].ay;
        p_h[i].vz += _dt*particles[i].az;
    }
}

static void reb_whfasthelio_keplerstep(const struct reb_simulation* const r, const double _dt){
    const int N_real = r->N-r->N_var;
    struct reb_particle* const p_h = r->ri_whfasthelio.p_h;
    const double M = r->particles[0].m*r->G;
#pragma omp parallel for
    for (unsigned int i=1;i<N_real;i++){
        kepler_step(r, p_h, M, i, _dt);
    }
    p_h[0].x += _dt*p_h[0].vx;
    p_h[0].y += _dt*p_h[0].vy;
    p_h[0].z += _dt*p_h[0].vz;
}

/***************************** 
 * Correctors                */

static void reb_whfasthelio_corrector_Z(struct reb_simulation* r, const double a, const double b){
    struct reb_simulation_integrator_whfasthelio* const ri_whfasthelio = &(r->ri_whfasthelio);
    struct reb_particle* restrict const particles = r->particles;
    const int N_real = r->N-r->N_var;
    reb_whfasthelio_keplerstep(r, a);
    to_inertial_pos(particles, ri_whfasthelio->p_h, N_real);
    reb_update_acceleration(r);
    reb_whfasthelio_interaction_step(r,-b);
    reb_whfasthelio_jump_step(r,-b);
    reb_whfasthelio_keplerstep(r, -2.*a);
    to_inertial_pos(particles, ri_whfasthelio->p_h, N_real);
    reb_update_acceleration(r);
    reb_whfasthelio_interaction_step(r,b);
    reb_whfasthelio_jump_step(r,b);
    reb_whfasthelio_keplerstep(r, a);
}


void reb_integrator_whfasthelio_part1(struct reb_simulation* const r){
    if (r->var_config_N){
        reb_exit("WHFastHELIO does currently not work with variational equations.");
    }
    struct reb_simulation_integrator_whfasthelio* const ri_whfasthelio = &(r->ri_whfasthelio);
    struct reb_particle* restrict const particles = r->particles;
    const int N_real = r->N - r->N_var;
    r->gravity_ignore_terms = 2;


    if (ri_whfasthelio->allocated_N != N_real){
        ri_whfasthelio->allocated_N = N_real;
        ri_whfasthelio->p_h = realloc(ri_whfasthelio->p_h,sizeof(struct reb_particle)*N_real);
        ri_whfasthelio->recalculate_heliocentric_this_timestep = 1;
    }

    if (ri_whfasthelio->safe_mode || ri_whfasthelio->recalculate_heliocentric_this_timestep == 1){
        if (ri_whfasthelio->is_synchronized==0){
            reb_integrator_whfasthelio_synchronize(r);
            if (ri_whfasthelio->recalculate_heliocentric_but_not_synchronized_warning==0){
                reb_warning(r,"Recalculating heliocentric coordinates but pos/vel were not synchronized before.");
                ri_whfasthelio->recalculate_heliocentric_but_not_synchronized_warning++;
            }
        }
        ri_whfasthelio->recalculate_heliocentric_this_timestep = 0;
        to_helio_posvel(particles, ri_whfasthelio->p_h, N_real);
    }

    if (ri_whfasthelio->is_synchronized==1){
        // First half DRIFT step
        if (ri_whfasthelio->corrector){
            reb_whfast_apply_corrector(r, 1.,ri_whfasthelio->corrector,reb_whfasthelio_corrector_Z);
        }
        reb_whfasthelio_keplerstep(r,r->dt/2.);
    }else{
        // Combined DRIFT step
        reb_whfasthelio_keplerstep(r,r->dt);
    }
    
    // For force calculation:
    if (r->force_is_velocity_dependent){
        to_inertial_posvel(particles, ri_whfasthelio->p_h, N_real);
    }else{
        to_inertial_pos(particles, ri_whfasthelio->p_h, N_real);
    }

    r->t+=r->dt/2.;
}

void reb_integrator_whfasthelio_synchronize(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfasthelio* const ri_whfasthelio = &(r->ri_whfasthelio);
    if (ri_whfasthelio->is_synchronized==0){
        const int N_real = r->N - r->N_var;
        struct reb_particle* sync_ph  = NULL;
        if (ri_whfasthelio->keep_unsynchronized){
            sync_ph = malloc(sizeof(struct reb_particle)*r->N);
            memcpy(sync_ph,r->ri_whfasthelio.p_h,r->N*sizeof(struct reb_particle));
        }
        struct reb_particle* restrict const particles = r->particles;
        reb_whfasthelio_keplerstep(r,r->dt/2.);
        if (ri_whfasthelio->corrector){
            reb_whfast_apply_corrector(r, -1.,ri_whfasthelio->corrector,reb_whfasthelio_corrector_Z);
        }
        to_inertial_posvel(particles, ri_whfasthelio->p_h, N_real);
        if (ri_whfasthelio->keep_unsynchronized){
            memcpy(r->ri_whfasthelio.p_h,sync_ph,r->N*sizeof(struct reb_particle));
            free(sync_ph);
        }else{
            ri_whfasthelio->is_synchronized=1;
        }
    }
}

void reb_integrator_whfasthelio_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfasthelio* const ri_whfasthelio = &(r->ri_whfasthelio);

    reb_whfasthelio_interaction_step(r,r->dt);
    reb_whfasthelio_jump_step(r,r->dt);
    
    ri_whfasthelio->is_synchronized=0;
    if (ri_whfasthelio->safe_mode){
        reb_integrator_whfasthelio_synchronize(r);
    }

    r->t+=r->dt/2.;
    r->dt_last_done = r->dt;
}
    
void reb_integrator_whfasthelio_reset(struct reb_simulation* const r){
    struct reb_simulation_integrator_whfasthelio* const ri_whfasthelio = &(r->ri_whfasthelio);
    ri_whfasthelio->allocated_N = 0;
    ri_whfasthelio->safe_mode = 1;
    ri_whfasthelio->recalculate_heliocentric_this_timestep = 0;
    ri_whfasthelio->recalculate_heliocentric_but_not_synchronized_warning = 0;
    ri_whfasthelio->is_synchronized = 1;
    if (ri_whfasthelio->p_h){
        free(ri_whfasthelio->p_h);
        ri_whfasthelio->p_h = NULL;
    }
}
