
/**
 * @file    integrator_whfasthelio.c
 * @brief   WHFASTHELIO integration scheme.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the WHFast integration scheme in
 *          Democratic Heliocentric Coordinates.
 *          It uses the alternative splitting WHDS proposed by Hernandez 
 *          and Dehnen (2017) which splits the Hamiltonian into
 *          three operators.
 *          The Kepler Solver is the same as in WHFast, 
 *          described in Rein & Tamayo 2015.
 * 
 * @section LICENSE
 * Copyright (c) 2016 Hanno Rein, Daniel Tamayo, Ari Silburt
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

/***************************** 
 * Operators                 */
static void reb_whfasthelio_jump_step(const struct reb_simulation* const r, double _dt){
    const int N_real = r->N-r->N_var;
    struct reb_particle* const p_h = r->ri_whfasthelio.p_h;
    const double m0 = r->particles[0].m;
    double px=0, py=0, pz=0;
    for(int i=1;i<N_real;i++){
        const double m = r->particles[i].m;
        px += m * p_h[i].vx / (m0+m);
        py += m * p_h[i].vy / (m0+m);
        pz += m * p_h[i].vz / (m0+m);
    }
    for(int i=1;i<N_real;i++){
        const double m = r->particles[i].m;
        p_h[i].x += _dt * (px - (m * p_h[i].vx / (m0+m)) );
        p_h[i].y += _dt * (py - (m * p_h[i].vy / (m0+m)) );
        p_h[i].z += _dt * (pz - (m * p_h[i].vz / (m0+m)) );
    }
}

static void reb_whfasthelio_interaction_step(const struct reb_simulation* const r, const double _dt){
    struct reb_particle* particles = r->particles;
    const int N_real = r->N-r->N_var;
    struct reb_particle* const p_h = r->ri_whfasthelio.p_h;
    const double m0 = r->particles[0].m;   
    for (unsigned int i=1;i<N_real;i++){
        const double m = r->particles[i].m;  
        p_h[i].vx += _dt*particles[i].ax*(m+m0)/m0;
        p_h[i].vy += _dt*particles[i].ay*(m+m0)/m0;
        p_h[i].vz += _dt*particles[i].az*(m+m0)/m0;
    }
}

static void reb_whfasthelio_keplerstep(const struct reb_simulation* const r, const double _dt){
    const int N_real = r->N-r->N_var;
    struct reb_particle* const p_h = r->ri_whfasthelio.p_h;
    const double m0 = r->particles[0].m;
#pragma omp parallel for
    for (unsigned int i=1;i<N_real;i++){
        kepler_step(r, p_h, r->G*(p_h[i].m + m0), i, _dt);
    }
    p_h[0].x += _dt*p_h[0].vx;
    p_h[0].y += _dt*p_h[0].vy;
    p_h[0].z += _dt*p_h[0].vz;
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
        reb_transformations_inertial_to_democratic_heliocentric_posvel(particles, ri_whfasthelio->p_h, N_real);
    }

    if (ri_whfasthelio->is_synchronized==1){
        // First half DRIFT step
        reb_whfasthelio_keplerstep(r,r->dt/2.);
    }else{
        // Combined DRIFT step
        reb_whfasthelio_keplerstep(r,r->dt);
    }
    
    reb_whfasthelio_jump_step(r,r->dt/2.);

    // For force calculation:
    if (r->force_is_velocity_dependent){
        reb_transformations_democratic_heliocentric_to_inertial_posvel(particles, ri_whfasthelio->p_h, N_real);
    }else{
        reb_transformations_democratic_heliocentric_to_inertial_pos(particles, ri_whfasthelio->p_h, N_real);
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
        reb_transformations_democratic_heliocentric_to_inertial_posvel(particles, ri_whfasthelio->p_h, N_real);
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
    reb_whfasthelio_jump_step(r,r->dt/2.);
    
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
