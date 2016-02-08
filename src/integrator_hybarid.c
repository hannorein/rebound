/**
 * @file 	integrator_hybarid.c
 * @brief 	WHFAST/IAS15 Hybrid integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements a hybrid integration scheme capable
 *  of handling close encounters, simple collisions, and
 *  planetesimal forces.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "output.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b

static void reb_integrator_hybarid_check_for_encounter(struct reb_simulation* r);
static void reb_integrator_hybarid_additional_forces_mini(struct reb_simulation* mini);


void reb_integrator_hybarid_part1(struct reb_simulation* r){
	const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    if (r->ri_hybarid.mini == NULL){
        r->ri_hybarid.mini = reb_create_simulation();
        r->ri_hybarid.mini->usleep = -1; // Disable visualiation
        r->ri_hybarid.mini->integrator = REB_INTEGRATOR_IAS15;
        r->ri_hybarid.mini->additional_forces = reb_integrator_hybarid_additional_forces_mini;
        r->ri_hybarid.mini->ri_hybarid.global = r;
        r->ri_hybarid.mini->testparticle_type = r->testparticle_type;
        r->ri_hybarid.mini->collision = r->collision;
        r->ri_hybarid.mini->collision_resolve = r->collision_resolve;
        r->ri_hybarid.mini->ri_ias15.epsilon = r->ri_ias15.epsilon;   //A.S. TEMPP
        //r->ri_hybarid.mini->ri_ias15.epsilon = 1e-8;  //speeds up ias and hybarid immensely
    }

    // Remove all particles from mini
    r->ri_hybarid.mini->t = r->t;
    r->ri_hybarid.mini->N = 0;
    r->ri_hybarid.mini->N_active = -1;
    r->ri_hybarid.mini_active = 0;
    r->ri_hybarid.global_index_from_mini_index_N = 0;
    
    //reset is_in_mini
    if (r->N>r->ri_hybarid.is_in_mini_Nmax){
        r->ri_hybarid.is_in_mini_Nmax = r->N;
        r->ri_hybarid.is_in_mini = realloc(r->ri_hybarid.is_in_mini,r->N*sizeof(int));
    }
    for(int i=0;i<r->N;i++)r->ri_hybarid.is_in_mini[i] = 0;

    // Add all massive particles
    for (int i=0; i<_N_active; i++){
        reb_add(r->ri_hybarid.mini, r->particles[i]);
        r->ri_hybarid.is_in_mini[i] = 1;
        if (r->ri_hybarid.global_index_from_mini_index_N>=r->ri_hybarid.encounter_index_Nmax){
            r->ri_hybarid.encounter_index_Nmax += 32;
            r->ri_hybarid.global_index_from_mini_index = realloc(r->ri_hybarid.global_index_from_mini_index,r->ri_hybarid.encounter_index_Nmax*sizeof(int));
        }
        r->ri_hybarid.global_index_from_mini_index[r->ri_hybarid.global_index_from_mini_index_N] = i;
        r->ri_hybarid.global_index_from_mini_index_N++;
    }
    r->ri_hybarid.mini->N_active = _N_active;

    reb_integrator_hybarid_check_for_encounter(r);
    
    //keep this after check_for_encounter - then no need to re-edit particles_prev
    if (r->testparticle_type && r->ri_hybarid.mini_active){
        if (r->N>r->ri_hybarid.particles_prev_Nmax){
            r->ri_hybarid.particles_prev_Nmax = r->N;
            r->ri_hybarid.particles_prev = realloc(r->ri_hybarid.particles_prev,r->N*sizeof(struct reb_particle));
        }
        memcpy(r->ri_hybarid.particles_prev, r->particles, sizeof(struct reb_particle)*r->N);
    }
    
    reb_integrator_whfast_part1(r);
}

void reb_integrator_hybarid_part2(struct reb_simulation* r){
    reb_integrator_whfast_part2(r);

    struct reb_simulation* mini = r->ri_hybarid.mini;
    if (r->ri_hybarid.mini_active){
        reb_integrate(mini,r->t);
        for (int i=0; i<mini->N; i++){
            r->particles[r->ri_hybarid.global_index_from_mini_index[i]] = mini->particles[i];
            r->particles[r->ri_hybarid.global_index_from_mini_index[i]].sim = r;
        }
        // Correct for energy jump in collision
        if (r->ri_hybarid.collision_in_timestep!=0){
            double Ei = r->ri_hybarid.energy_before_collision_in_timestep;
            double Ef = reb_tools_energy(r);
            r->ri_hybarid.dE_offset += Ei - Ef;
            r->ri_hybarid.collision_in_timestep = 0;
        }
    }
}
	
void reb_integrator_hybarid_synchronize(struct reb_simulation* r){
	// Do nothing.
    reb_integrator_whfast_synchronize(r);
}

void reb_integrator_hybarid_reset(struct reb_simulation* r){
    r->ri_hybarid.timestep_too_large_warning = 0.;
	// TODO: Implement rest
    reb_integrator_whfast_reset(r);
}

static void reb_integrator_hybarid_check_for_encounter(struct reb_simulation* global){
    struct reb_simulation* mini = global->ri_hybarid.mini;
	const int _N_active = ((global->N_active==-1)?global->N:global->N_active) - global->N_var;
    struct reb_particle* global_particles = global->particles;
    struct reb_particle p0 = global_particles[0];
    double switch_ratio = global->ri_hybarid.switch_ratio;
    double switch_ratio2 = switch_ratio*switch_ratio;
    double min_dt_enc2 = INFINITY;
    for (int i=0; i<_N_active; i++){
        struct reb_particle pi = global_particles[i];
        double radius_check2 = global->ri_hybarid.CE_radius*pi.r;
        const double dxi = p0.x - pi.x;
        const double dyi = p0.y - pi.y;
        const double dzi = p0.z - pi.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double mi = pi.m/(p0.m*3.);
        double rhi = pow(mi*mi*r0i2*r0i2*r0i2,1./6.);
        for(int j=i+1;j<global->N;j++){
            struct reb_particle pj = global_particles[j];
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double mj = pj.m/(p0.m*3.);
            double rhj = pow(mj*mj*r0j2*r0j2*r0j2,1./6.);
            const double rh_sum = rhi+rhj;
            const double rh_sum2 = rh_sum*rh_sum;
            
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;

            if(rij2 < switch_ratio2*rh_sum2 || rij2 < radius_check2){
                global->ri_hybarid.mini_active = 1;
                if (j>=_N_active && global->ri_hybarid.is_in_mini[j]==0){//make sure not already added
                    // Monitor hillradius/relative velocity
                    const double dvx = pi.vx - pj.vx;
                    const double dvy = pi.vy - pj.vy;
                    const double dvz = pi.vz - pj.vz;
                    const double vij2 = dvx*dvx + dvy*dvy + dvz*dvz;
                    const double dt_enc2 = switch_ratio2*rh_sum2/vij2;
                    min_dt_enc2 = MIN(min_dt_enc2,dt_enc2);
                    // Add particle to mini simulation
                    reb_add(mini,pj);
                    global->ri_hybarid.is_in_mini[j] = 1;
                    if (global->ri_hybarid.global_index_from_mini_index_N>=global->ri_hybarid.encounter_index_Nmax){
                        while(global->ri_hybarid.global_index_from_mini_index_N>=global->ri_hybarid.encounter_index_Nmax) global->ri_hybarid.encounter_index_Nmax += 32;
                        global->ri_hybarid.global_index_from_mini_index = realloc(global->ri_hybarid.global_index_from_mini_index,global->ri_hybarid.encounter_index_Nmax*sizeof(int));
                    }
                    global->ri_hybarid.global_index_from_mini_index[global->ri_hybarid.global_index_from_mini_index_N] = j;
                    global->ri_hybarid.global_index_from_mini_index_N++;
                }
            }
        }
    }
    if (global->ri_hybarid.timestep_too_large_warning==0 && min_dt_enc2 < 16.*global->dt*global->dt){
        global->ri_hybarid.timestep_too_large_warning = 1;
        reb_warning("The timestep is likely too large. Close encounters might be missed. Decrease the timestep or increase the switching radius. This warning will appear only once.");
    }
}

static void reb_integrator_hybarid_additional_forces_mini(struct reb_simulation* mini){
    if (mini->testparticle_type){
        struct reb_simulation* global = mini->ri_hybarid.global;
        struct reb_particle* global_particles = global->particles;
        struct reb_particle* mini_particles = mini->particles;
        struct reb_particle* global_prev = global->ri_hybarid.particles_prev;
        const double t_prev = global->t - global->dt;
        double timefac = (mini->t - t_prev)/global->dt;
        // TODO: See if the following is good enough and if so why
        /// timefac = 0.5;
        const int N_active = global->N_active;
        const double G = global->G;
        for(int j=N_active;j<global->N;j++){            //planetesimals in global
            if(global->ri_hybarid.is_in_mini[j]==0){
                for(int i=0;i<mini->N_active;i++){      //massive bodies in mini
                    struct reb_particle* body = &(mini_particles[i]);
                    const double ix = (1.-timefac)*global_prev[j].x + timefac*global_particles[j].x; //interpolated values
                    const double iy = (1.-timefac)*global_prev[j].y + timefac*global_particles[j].y;
                    const double iz = (1.-timefac)*global_prev[j].z + timefac*global_particles[j].z;
                    const double mp = global_particles[j].m;
                    const double ddx = body->x - ix;
                    const double ddy = body->y - iy;
                    const double ddz = body->z - iz;
                    
                    const double rijinv2 = 1.0/(ddx*ddx + ddy*ddy + ddz*ddz);
                    const double ac = -G*mp*rijinv2*sqrt(rijinv2);
                    
                    body->ax += ac*ddx;     //perturbation on planets due to planetesimals.
                    body->ay += ac*ddy;
                    body->az += ac*ddz;
                }
            }
        }
    }
}

