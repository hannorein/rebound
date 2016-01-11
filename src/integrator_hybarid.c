/**
 * @file 	integrator.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
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
#include "integrator_ias15.h"
#include "integrator_whfast.h"

static void reb_integrator_hybarid_check_for_encounter(struct reb_simulation* r);


void reb_integrator_hybarid_additional_forces_mini(struct reb_simulation* mini){
    const double G = mini->G;
    const int N = mini->N;
    const int N_active = mini->N_active;
    struct reb_simulation* global = mini->ri_hybarid.global;
    struct reb_particle* particles_mini = mini->particles;
    struct reb_particle* particles_global = global->particles;
    struct reb_particle* particles_global_prev = global->ri_hybarid.particles_prev;
    double planetesimal_mass = 1e-8;
    
    const double Gm1 = G*planetesimal_mass;
    for(int i=0;i<N_active;i++){
        struct reb_particle* body = &(particles_mini[i]);
        for(int j=N_active;j<N;j++){//add planetesimal forces to massive bodies
            struct reb_particle p = particles_mini[j];
            
            const double dx = body->x - p.x;
            const double dy = body->y - p.y;
            const double dz = body->z - p.z;
            
            const double rijinv2 = 1.0/(dx*dx + dy*dy + dz*dz);
            const double ac = -Gm1*rijinv2*sqrt(rijinv2);
            
            body->ax += ac*dx;      //perturbation on planets due to planetesimals.
            body->ay += ac*dy;
            body->az += ac*dz;
        }
    }
    
    //forces from global into mini
    double t_prev = global->t - global->dt;
    const double timefac = (mini->t - t_prev)/global->dt;
    for(int i=N_active;i<global->N;i++){    //planetesimals
        if(global->ri_hybarid.is_in_mini[i]==0){             //find planetesimals which is part of global but not mini
            const double ix = timefac*particles_global[i].x - (1.-timefac)*particles_global_prev[i].x; //interpolated values
            const double iy = timefac*particles_global[i].y - (1.-timefac)*particles_global_prev[i].y;
            const double iz = timefac*particles_global[i].z - (1.-timefac)*particles_global_prev[i].z;
            for(int j=0;j<N_active;j++){//massive bodies
                struct reb_particle* body = &(particles_mini[j]);
                const double ddx = body->x - ix;
                const double ddy = body->y - iy;
                const double ddz = body->z - iz;
                
                const double rijinv2 = 1.0/(ddx*ddx + ddy*ddy + ddz*ddz);
                const double ac = -Gm1*rijinv2*sqrt(rijinv2);
                
                body->ax += ac*ddx;     //perturbation on planets due to planetesimals.
                body->ay += ac*ddy;
                body->az += ac*ddz;
            }
        }
    }

}

void reb_integrator_hybarid_part1(struct reb_simulation* r){
    if (r->ri_hybarid.mini == NULL){
        r->ri_hybarid.mini = reb_create_simulation();
        r->ri_hybarid.mini->integrator = REB_INTEGRATOR_IAS15;
        r->ri_hybarid.mini->additional_forces = reb_integrator_hybarid_additional_forces_mini;
        r->ri_hybarid.mini->ri_hybarid.global = r;
    }

    if (r->N>=r->ri_hybarid.particles_prev_Nmax){
        r->ri_hybarid.particles_prev_Nmax += 32;
        r->ri_hybarid.particles_prev = realloc(r->ri_hybarid.particles_prev,r->ri_hybarid.particles_prev_Nmax*sizeof(struct reb_particle));
    }
    memcpy(r->ri_hybarid.particles_prev, r->particles, sizeof(struct reb_particle)*r->N); 

    // Remove all particles from mini
    r->ri_hybarid.mini->t = r->t;
    r->ri_hybarid.mini->N = 0;
    r->ri_hybarid.mini->N_active = -1;
    r->ri_hybarid.mini_active = 0;
    r->ri_hybarid.encounter_index_N = 0;
    
    
    if (r->N>=r->ri_hybarid.is_in_mini_Nmax){
        r->ri_hybarid.is_in_mini_Nmax += 32;
        r->ri_hybarid.is_in_mini = realloc(r->ri_hybarid.is_in_mini,r->ri_hybarid.is_in_mini_Nmax*sizeof(int));
    }

    // Add all massive particles
    for (int i=0; i<r->N_active; i++){
        reb_add(r->ri_hybarid.mini, r->particles[i]);
        r->ri_hybarid.is_in_mini[i] = 1;
        if (r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax){
            r->ri_hybarid.encounter_index_Nmax += 32;
            r->ri_hybarid.encounter_index = realloc(r->ri_hybarid.encounter_index,r->ri_hybarid.encounter_index_Nmax*sizeof(int));
        }
        r->ri_hybarid.encounter_index[r->ri_hybarid.encounter_index_N] = i;
        r->ri_hybarid.encounter_index_N++;
    }

    reb_integrator_hybarid_check_for_encounter(r);

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


    // Check for encounters
    //   if new encounters, copy from global to mini
    //
    // Copy positions from mini to global
    //
    // Store the pos to prev_pos
    //
}
	
void reb_integrator_hybarid_synchronize(struct reb_simulation* r){
	// Do nothing.
    reb_integrator_whfast_synchronize(r);
}

void reb_integrator_hybarid_reset(struct reb_simulation* r){
	// Do nothing.
    reb_integrator_whfast_reset(r);
}


//collect the id/array number of all planetesimals involved in a close encounter
static void reb_integrator_hybarid_check_for_encounter(struct reb_simulation* r){
    struct reb_simulation* mini = r->ri_hybarid.mini;
    const int N = r->N;
    const int N_active = r->N_active;
    struct reb_particle* global = r->particles;
    struct reb_particle p0 = global[0];
    for (int i=0; i<N_active; i++){
        struct reb_particle pi = global[i];
        const double dxi = p0.x - pi.x;
        const double dyi = p0.y - pi.y;
        const double dzi = p0.z - pi.z;
        const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double rhi = r0i2*pow(pi.m/(p0.m*3.),2./3.);
        
        for (int j=i+1; j<N; j++){
            struct reb_particle pj = global[j];
            double HSR = r->ri_hybarid.switch_ratio;
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*pow(pj.m/(p0.m*3.),2./3.);
            
            const double dx = pi.x - pj.x;
            const double dy = pi.y - pj.y;
            const double dz = pi.z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double ratio = rij2/(rhi+rhj);    //(p-p distance/Hill radii)^2
            r->ri_hybarid.is_in_mini[j] = 0;
            
            if(ratio < HSR){
                printf("\nencounter found\n");
                r->ri_hybarid.mini_active = 1;
                if (j>=r->N_active){
                    reb_add(mini,pj);
                    r->ri_hybarid.is_in_mini[j] = 1;
                    if (r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax){
                        r->ri_hybarid.encounter_index_Nmax += 32;
                        r->ri_hybarid.encounter_index = realloc(r->ri_hybarid.encounter_index,r->ri_hybarid.encounter_index_Nmax*sizeof(int));
                    }
                    r->ri_hybarid.encounter_index[r->ri_hybarid.encounter_index_N] = j;
                    r->ri_hybarid.encounter_index_N++;
                }
            }
        }
    }
}
