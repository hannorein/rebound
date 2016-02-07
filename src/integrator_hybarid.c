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
    r->ri_hybarid.encounter_index_N = 0;
    
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
        if (r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax){
            r->ri_hybarid.encounter_index_Nmax += 32;
            r->ri_hybarid.encounter_index = realloc(r->ri_hybarid.encounter_index,r->ri_hybarid.encounter_index_Nmax*sizeof(int));
        }
        r->ri_hybarid.encounter_index[r->ri_hybarid.encounter_index_N] = i;
        r->ri_hybarid.encounter_index_N++;
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
            r->particles[r->ri_hybarid.encounter_index[i]] = mini->particles[i];
            r->particles[r->ri_hybarid.encounter_index[i]].sim = r;
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

static void reb_integrator_hybarid_check_for_encounter(struct reb_simulation* r){
    struct reb_simulation* mini = r->ri_hybarid.mini;
	const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    struct reb_particle* global = r->particles;
    struct reb_particle p0 = global[0];
    double HSR2 = r->ri_hybarid.switch_ratio*r->ri_hybarid.switch_ratio;
    double minr2v2 = INFINITY;
    for (int i=0; i<_N_active; i++){
        struct reb_particle* pi = &(global[i]);
        double radius_check2 = r->ri_hybarid.CE_radius*r->ri_hybarid.CE_radius*pi->r*pi->r;
        const double dxi = p0.x - pi->x;
        const double dyi = p0.y - pi->y;
        const double dzi = p0.z - pi->z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        double rhi = r0i2*pow(pi->m/(p0.m*3.),2./3.);
        for(int j=i+1;j<r->N;j++){
            struct reb_particle pj = global[j];
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*pow(pj.m/(p0.m*3.),2./3.);
            
            const double dx = pi->x - pj.x;
            const double dy = pi->y - pj.y;
            const double dz = pi->z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double ratio = rij2/(rhi+rhj);    //(p-p distance/Hill radii)^2

            if(ratio < HSR2 || rij2 < radius_check2){
                r->ri_hybarid.mini_active = 1;
                if (j>=_N_active && r->ri_hybarid.is_in_mini[j]==0){//make sure not already added
                    // Monitor hillradius/relative velocity
                    const double dvx = pi->vx - pj.vx;
                    const double dvy = pi->vy - pj.vy;
                    const double dvz = pi->vz - pj.vz;
                    const double vij2 = dvx*dvx + dvy*dvy + dvz*dvz;
                    const double r2v2 = HSR2*(rhi+rhj)*(rhi+rhj)/vij2;
                    minr2v2 = MIN(minr2v2,r2v2);
                    // Add particle to mini simulation
                    reb_add(mini,pj);
                    r->ri_hybarid.is_in_mini[j] = 1;
                    if (r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax){
                        while(r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax) r->ri_hybarid.encounter_index_Nmax += 32;
                        r->ri_hybarid.encounter_index = realloc(r->ri_hybarid.encounter_index,r->ri_hybarid.encounter_index_Nmax*sizeof(int));
                    }
                    r->ri_hybarid.encounter_index[r->ri_hybarid.encounter_index_N] = j;
                    r->ri_hybarid.encounter_index_N++;
                }
            }
        }
    }
    if (r->ri_hybarid.timestep_too_large_warning==0 && minr2v2<16.*r->dt*r->dt){
        r->ri_hybarid.timestep_too_large_warning = 1;
        reb_warning("The timestep appears to be too large. Close encounters might be missed. Decrease the timestep or increase the switching radius. This warning appear only once.");
    }
}

static void reb_integrator_hybarid_additional_forces_mini(struct reb_simulation* mini){
    if (mini->testparticle_type){
        struct reb_simulation* r = mini->ri_hybarid.global;
        struct reb_particle* global = r->particles;
        struct reb_particle* mini_particles = mini->particles;
        struct reb_particle* global_prev = r->ri_hybarid.particles_prev;
        const double t_prev = r->t - r->dt;
        const double timefac = (mini->t - t_prev)/r->dt;
        const int N_active = r->N_active;
        const double G = r->G;
        for(int i=N_active;i<r->N;i++){    //planetesimals
            if(r->ri_hybarid.is_in_mini[i]==0){
                const double ix = (1.-timefac)*global_prev[i].x + timefac*global[i].x; //interpolated values
                const double iy = (1.-timefac)*global_prev[i].y + timefac*global[i].y;
                const double iz = (1.-timefac)*global_prev[i].z + timefac*global[i].z;
                const double mp = global[i].m;
                for(int j=0;j<N_active;j++){//massive bodies
                    struct reb_particle* body = &(mini_particles[j]);
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

