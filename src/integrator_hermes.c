/**
 * @file 	integrator_hermes.c
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

static void reb_integrator_hermes_check_for_encounter(struct reb_simulation* r);
static void reb_integrator_hermes_additional_forces_mini(struct reb_simulation* mini);
static void calc_forces_on_planets(const struct reb_simulation* r, double* a);

void reb_integrator_hermes_part1(struct reb_simulation* r){
	const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    struct reb_simulation* mini = r->ri_hermes.mini;
    if (mini == NULL){
        mini = reb_create_simulation();
        r->ri_hermes.mini = mini;
        mini->usleep = -1; // Disable visualiation
        mini->integrator = REB_INTEGRATOR_IAS15;
        mini->gravity = REB_GRAVITY_BASIC;
        mini->dt = r->dt;
        mini->additional_forces = reb_integrator_hermes_additional_forces_mini;
        mini->ri_hermes.global = r;    //set to != 0 so that collision.c knows to remove from both
    }
    mini->testparticle_type = r->testparticle_type;
    mini->collision = r->collision;
    mini->collision_resolve = r->collision_resolve;
    mini->collision_resolve_keep_sorted = r->collision_resolve_keep_sorted;
    mini->track_energy_offset = r->track_energy_offset;

    // Remove all particles from mini
    mini->t = r->t;
    int mini_previously_active = r->ri_hermes.mini_active;
    mini->N = 0;
    mini->energy_offset = 0.;
    r->ri_hermes.mini_active = 0;
    r->ri_hermes.global_index_from_mini_index_N = 0;
    r->ri_hermes.collision_this_global_dt = 0;
    
    if (_N_active>r->ri_hermes.a_Nmax){
        r->ri_hermes.a_i = realloc(r->ri_hermes.a_i,sizeof(double)*3*_N_active);
        r->ri_hermes.a_f = realloc(r->ri_hermes.a_f,sizeof(double)*3*_N_active);
        r->ri_hermes.a_Nmax = _N_active;
    }
    
    //reset is_in_mini
    if (r->N>r->ri_hermes.is_in_mini_Nmax){
        r->ri_hermes.is_in_mini_Nmax = r->N;
        r->ri_hermes.is_in_mini = realloc(r->ri_hermes.is_in_mini,r->N*sizeof(int));
    }
    for(int i=_N_active;i<r->N;i++)r->ri_hermes.is_in_mini[i] = 0;
    
    // Add all massive particles
    for (int i=0; i<_N_active; i++){
        reb_add(r->ri_hermes.mini, r->particles[i]);
        r->ri_hermes.is_in_mini[i] = 1;
        if (r->ri_hermes.global_index_from_mini_index_N>=r->ri_hermes.global_index_from_mini_index_Nmax){
            r->ri_hermes.global_index_from_mini_index_Nmax += 32;
            r->ri_hermes.global_index_from_mini_index = realloc(r->ri_hermes.global_index_from_mini_index,r->ri_hermes.global_index_from_mini_index_Nmax*sizeof(int));
        }
        r->ri_hermes.global_index_from_mini_index[r->ri_hermes.global_index_from_mini_index_N] = i;
        r->ri_hermes.global_index_from_mini_index_N++;
    }
    r->ri_hermes.mini->N_active = _N_active;

    reb_integrator_hermes_check_for_encounter(r);
        
    if (r->N != r->ri_hermes.mini->N || mini_previously_active==0) {
        reb_integrator_ias15_clear(r->ri_hermes.mini);
    }
    
    calc_forces_on_planets(r, r->ri_hermes.a_i);
    
    if(r->ri_hermes.mini_active && r->track_energy_offset){
        r->ri_hermes.energy_before_timestep = reb_tools_energy(r);
    }
    
    reb_integrator_whfast_part1(r);
}


void reb_integrator_hermes_part2(struct reb_simulation* r){
    reb_integrator_whfast_part2(r);
    
    calc_forces_on_planets(r, r->ri_hermes.a_f);
    
    struct reb_simulation* mini = r->ri_hermes.mini;
    r->ri_hermes.steps++;
    if (r->ri_hermes.mini_active){
        r->ri_hermes.steps_miniactive++;
        r->ri_hermes.steps_miniN += mini->N;
        reb_integrate(mini,r->t);

        for (int i=0; i<mini->N; i++){
            r->particles[r->ri_hermes.global_index_from_mini_index[i]] = mini->particles[i];
            r->particles[r->ri_hermes.global_index_from_mini_index[i]].sim = r;    
        }
        
        // Correct for energy jump in collision
        if(r->ri_hermes.collision_this_global_dt && r->track_energy_offset){
            double Ef = reb_tools_energy(r);
            r->energy_offset += r->ri_hermes.energy_before_timestep - Ef;
        }
    }
}
	
void reb_integrator_hermes_synchronize(struct reb_simulation* r){
	// Do nothing.
    reb_integrator_whfast_synchronize(r);
}

void reb_integrator_hermes_reset(struct reb_simulation* r){
    //r->ri_hermes.timestep_too_large_warning = 0.; //Don't think we want to reset the warning.
    r->ri_hermes.steps = 0;
    r->ri_hermes.steps_miniactive = 0;
    r->ri_hermes.steps_miniN = 0;
    
    reb_integrator_whfast_reset(r);

    if (r->ri_hermes.mini){
        reb_free_simulation(r->ri_hermes.mini);
        r->ri_hermes.mini = NULL;
    }
    if(r->ri_hermes.global_index_from_mini_index){
        free(r->ri_hermes.global_index_from_mini_index);
        r->ri_hermes.global_index_from_mini_index = NULL;
        r->ri_hermes.global_index_from_mini_index_Nmax = 0;
    }
    if(r->ri_hermes.is_in_mini){
        free(r->ri_hermes.is_in_mini);
        r->ri_hermes.is_in_mini = NULL;
        r->ri_hermes.is_in_mini_Nmax = 0;
    }
    if(r->ri_hermes.a_i){
        free(r->ri_hermes.a_i);
    }
    if(r->ri_hermes.a_f){
        free(r->ri_hermes.a_f);
    }
    r->ri_hermes.a_Nmax = 0;
}

static void reb_integrator_hermes_check_for_encounter(struct reb_simulation* global){
    struct reb_simulation* mini = global->ri_hermes.mini;
	const int _N_active = ((global->N_active==-1)?global->N:global->N_active) - global->N_var;
    struct reb_particle* global_particles = global->particles;
    struct reb_particle p0 = global_particles[0];
    double hill_switch_factor = global->ri_hermes.hill_switch_factor;
    double hill_switch_factor2 = hill_switch_factor*hill_switch_factor;
    double min_dt_enc2 = INFINITY;
    for (int i=0; i<_N_active; i++){
        struct reb_particle pi = global_particles[i];
        double radius_check = global->ri_hermes.radius_switch_factor*pi.r;
        double radius_check2 = radius_check*radius_check;
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

            if(rij2 < hill_switch_factor2*rh_sum2 || rij2 < radius_check2){
                global->ri_hermes.mini_active = 1;
                // Monitor hill radius/relative velocity
                const double dvx = pi.vx - pj.vx;
                const double dvy = pi.vy - pj.vy;
                const double dvz = pi.vz - pj.vz;
                const double vij2 = dvx*dvx + dvy*dvy + dvz*dvz;
                const double dt_enc2 = hill_switch_factor2*rh_sum2/vij2;
                min_dt_enc2 = MIN(min_dt_enc2,dt_enc2);
                if (j>=_N_active && global->ri_hermes.is_in_mini[j]==0){//make sure not already added
                    // Add particle to mini simulation
                    reb_add(mini,pj);
                    global->ri_hermes.is_in_mini[j] = 1;
                    if (global->ri_hermes.global_index_from_mini_index_N>=global->ri_hermes.global_index_from_mini_index_Nmax){
                        while(global->ri_hermes.global_index_from_mini_index_N>=global->ri_hermes.global_index_from_mini_index_Nmax) global->ri_hermes.global_index_from_mini_index_Nmax += 32;
                        global->ri_hermes.global_index_from_mini_index = realloc(global->ri_hermes.global_index_from_mini_index,global->ri_hermes.global_index_from_mini_index_Nmax*sizeof(int));
                    }
                    global->ri_hermes.global_index_from_mini_index[global->ri_hermes.global_index_from_mini_index_N] = j;
                    global->ri_hermes.global_index_from_mini_index_N++;
                }
            }
        }
    }
    if (global->ri_hermes.timestep_too_large_warning==0 && min_dt_enc2 < 16.*global->dt*global->dt){
        global->ri_hermes.timestep_too_large_warning = 1;
        reb_warning("The timestep is likely too large. Close encounters might be missed. Decrease the timestep or increase the switching radius. This warning will appear only once.");
    }
}

static void calc_forces_on_planets(const struct reb_simulation* r, double* a){
    int* is_in_mini = r->ri_hermes.is_in_mini;
    double G = r->G;
    const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    for (int i = 0; i<_N_active; i++){
        struct reb_particle pm = r->particles[i];
        double ax = 0.;
        double ay = 0.;
        double az = 0.;
        for (int j = _N_active; j<r->N; j++){
            if (is_in_mini[j] == 0){
                struct reb_particle ps = r->particles[j];
                double dx = ps.x - pm.x;
                double dy = ps.y - pm.y;
                double dz = ps.z - pm.z;
                double d = sqrt(dx*dx + dy*dy + dz*dz);
                ax += ps.m * dx * G/(d*d*d);
                ay += ps.m * dy * G/(d*d*d);
                az += ps.m * dz * G/(d*d*d);
            }
        }
        a[i*3+0] = ax;
        a[i*3+1] = ay;
        a[i*3+2] = az;
    }
}

// This is the current algorithm, interpolating forces
static void reb_integrator_hermes_additional_forces_mini(struct reb_simulation* mini){
    if (mini->testparticle_type){
        struct reb_simulation* global = mini->ri_hermes.global;
        struct reb_particle* mini_particles = mini->particles;
        const double t_prev = global->t - global->dt;
        double timefac = (mini->t - t_prev)/global->dt;
        
        double* a_i = global->ri_hermes.a_i;
        double* a_f = global->ri_hermes.a_f;
        // TODO: See if the following is good enough and if so why
        // timefac = 0.5;
#pragma omp parallel for schedule(guided)
        for(int i=0;i<mini->N_active;i++){              //massive bodies in mini
            double ax0 = a_i[i*3+0];
            double ay0 = a_i[i*3+1];
            double az0 = a_i[i*3+2];
            double ax1 = a_f[i*3+0];
            double ay1 = a_f[i*3+1];
            double az1 = a_f[i*3+2];
            
            mini_particles[i].ax += ax0*(1.-timefac) + ax1*timefac;
            mini_particles[i].ay += ay0*(1.-timefac) + ay1*timefac;
            mini_particles[i].az += az0*(1.-timefac) + az1*timefac;
        }
    }
}
