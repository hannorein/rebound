/**
 * @file    integrator_hermes.c
 * @brief   HERMES. A WHFAST/IAS15 hybrid integration scheme.
 * @author  Ari Silburt <silburt@astro.utoronto.ca>
 * @details This file implements a hybrid integration scheme capable
 *  of handling close encounters, simple collisions, and
 *  planetesimal forces. Details are describe in Silburt et al (in prep).
 * 
 * @section LICENSE
 * Copyright (c) 2016 Ari Silburt 
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
#include "integrator_whfasthelio.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

static void reb_integrator_hermes_check_for_encounter(struct reb_simulation* r);
static void reb_integrator_hermes_additional_forces_mini(struct reb_simulation* mini);
static void reb_integrator_hermes_apply_forces(const struct reb_simulation* r, double* a);
static void reb_integrator_hermes_get_ae(struct reb_simulation* r, int index, double* a, double* e);
static void reb_integrator_hermes_autocalc_HSF(struct reb_simulation* r, double* min_dt_enc2, int i, int j);
static void reb_integrator_hermes_autocalc_HSF_case_total_overlap(double mu, double ri_min, double ri_max, double ai, double ei, double* vphi_max_i, double* vr_max_i, double aj, double ej, double rj_min, double rj_max, double* vphi_max_j, double* vr_max_j);
static void reb_integrator_hermes_autocalc_HSF_case_partial_overlap(double mu, double ri_min, double ri_max, double ai, double ei, double* vphi_max_i, double* vr_max_i, double aj, double ej, double rj_min, double rj_max, double* vphi_max_j, double* vr_max_j);

void reb_integrator_hermes_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
    const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    struct reb_simulation* mini = r->ri_hermes.mini;
    if (mini == NULL){
        mini = reb_create_simulation();
        r->ri_hermes.mini = mini;
        mini->visualization = REB_VISUALIZATION_NONE; // Disable visualiation
        mini->integrator = REB_INTEGRATOR_IAS15;
        mini->gravity = REB_GRAVITY_BASIC;
        mini->dt = r->dt;
        mini->additional_forces = reb_integrator_hermes_additional_forces_mini;
        mini->G = r->G;
        mini->softening = r->softening;
        if(r->collision_resolve_keep_sorted ==0) reb_warning(r,"When using HERMES, the user must set r->collision_resolve_keep_sorted = 1, or else it is likely that the wrong particle will be removed from the simulation during a collision/ejection, leading to energy jumps and other unpredictable behaviour. This warning will only appear once.\n");
    }
    mini->ri_hermes.global = r;    //set to != 0 so that collision.c knows to remove from both
    mini->testparticle_type = r->testparticle_type;
    mini->collision = r->collision;
    mini->collision_resolve = r->collision_resolve;
    mini->collision_resolve_keep_sorted = r->collision_resolve_keep_sorted;
    mini->track_energy_offset = r->track_energy_offset;
    mini->force_is_velocity_dependent = r->force_is_velocity_dependent;
    mini->post_timestep_modifications = r->post_timestep_modifications;

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
            while(r->ri_hermes.global_index_from_mini_index_N>=r->ri_hermes.global_index_from_mini_index_Nmax) r->ri_hermes.global_index_from_mini_index_Nmax += 32;
            r->ri_hermes.global_index_from_mini_index = realloc(r->ri_hermes.global_index_from_mini_index,r->ri_hermes.global_index_from_mini_index_Nmax*sizeof(int));
        }
        r->ri_hermes.global_index_from_mini_index[r->ri_hermes.global_index_from_mini_index_N] = i;
        r->ri_hermes.global_index_from_mini_index_N++;
    }
    r->ri_hermes.mini->N_active = _N_active;

    // Check for Close Encounters and Determine HSF
    r->ri_hermes.current_hill_switch_factor = r->ri_hermes.hill_switch_factor;
    reb_integrator_hermes_check_for_encounter(r);
        
    if (r->N != r->ri_hermes.mini->N || mini_previously_active==0) {
        reb_integrator_ias15_clear(r->ri_hermes.mini);
    }
    
    reb_integrator_hermes_apply_forces(r, r->ri_hermes.a_i);
    
    //reb_integrator_whfast_part1(r);
    reb_integrator_whfasthelio_part1(r);
}


void reb_integrator_hermes_part2(struct reb_simulation* r){
    //reb_integrator_whfast_part2(r);
    reb_integrator_whfasthelio_part2(r);
    
    reb_integrator_hermes_apply_forces(r, r->ri_hermes.a_f);
    
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
        if(r->track_energy_offset)r->energy_offset += r->ri_hermes.mini->energy_offset;

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
    double solar_check = global->ri_hermes.solar_switch_factor*p0.r;
    double solar_check2 = solar_check*solar_check;
    double current_hill_switch_factor = global->ri_hermes.current_hill_switch_factor;
    double hill_switch_factor2 = current_hill_switch_factor*current_hill_switch_factor;
    double min_dt_enc2 = INFINITY;
    double min_dt_enc2_autoHSF = INFINITY;
    for (int i=0; i<_N_active; i++){
        struct reb_particle pi = global_particles[i];
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

            if((rij2 < hill_switch_factor2*rh_sum2 && i>0) || (rij2 < solar_check2 && i==0)){
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
            } else if(global->ri_hermes.adaptive_hill_switch_factor) reb_integrator_hermes_autocalc_HSF(global, &min_dt_enc2_autoHSF, i, j); //find autoHSF
        }
    }
    if(global->ri_hermes.adaptive_hill_switch_factor){      //Calc optimal HSF value from min_dt_enc2_autoHSF value found in for loop
        if(min_dt_enc2_autoHSF < INFINITY){
            double dt2 = 16.*global->dt*global->dt;         //Factor of 4:hill sphere > 4*length scales for wiggle room
            double HSF_new = sqrt(dt2/min_dt_enc2_autoHSF);
            double base = 1.25;
            double exp = ceilf(log10(HSF_new)/log10(base)); //round HSF up to nearest multiple of 1.25
            HSF_new = pow(base,exp);
            global->ri_hermes.current_hill_switch_factor = MAX(global->ri_hermes.current_hill_switch_factor, HSF_new); // Increase HSF if needed
        }
    }
    
    if (global->ri_hermes.adaptive_hill_switch_factor==0 && global->ri_hermes.timestep_too_large_warning==0 && min_dt_enc2 < 16.*global->dt*global->dt){
        global->ri_hermes.timestep_too_large_warning = 1;
        reb_warning(global,"The timestep is likely too large. Close encounters might be missed. Decrease the timestep or increase the switching radius. This warning will appear only once.");
    }
}

//get min encounter time between overlapping orbits
//Using subscripts i and j for related variables, corresponding to bodies 'i' and 'j'
static void reb_integrator_hermes_autocalc_HSF(struct reb_simulation* r, double* min_dt_enc2, int i, int j){
    const double m0 = r->particles[0].m;
    const double mu = r->G*m0;
    double ei, ai, ej, aj;
    reb_integrator_hermes_get_ae(r, i, &ai, &ei);
    reb_integrator_hermes_get_ae(r, j, &aj, &ej);
    double ri_min = ai*(1.-ei);
    double ri_max = ai*(1.+ei);
    double rj_min = aj*(1.-ej);
    double rj_max = aj*(1.+ej);
    double vphi_max_i=0., vr_max_i=0.;            //max phi_hat velocity, max r_hat velocity for body i
    double vphi_max_j=0., vr_max_j=0.;            //max phi_hat velocity, max r_hat velocity for body j
    if((ri_max>rj_max)&&(ri_min<rj_min)){         //CASE1: body i totally overlaps body j
        reb_integrator_hermes_autocalc_HSF_case_total_overlap(mu, ri_min, ri_max, ai, ei, &vphi_max_i, &vr_max_i, aj, ej, rj_min, rj_max, &vphi_max_j, &vr_max_j);
    } else if((rj_max>ri_max)&&(rj_min<ri_min)){  //CASE2: body j totally overlaps body i
        reb_integrator_hermes_autocalc_HSF_case_total_overlap(mu, rj_min, rj_max, aj, ej, &vphi_max_j, &vr_max_j, ai, ei, ri_min, ri_max, &vphi_max_i, &vr_max_i);
    } else if((rj_max>ri_max)&&(rj_min<ri_max)){  //CASE3: partial overlap (i=inner body), boundaries: inner=rj_min, outer=ri_max
        reb_integrator_hermes_autocalc_HSF_case_partial_overlap(mu, ri_min, ri_max, ai, ei, &vphi_max_i, &vr_max_i, aj, ej, rj_min, rj_max, &vphi_max_j, &vr_max_j);
    } else if((ri_max>rj_max)&&(ri_min<rj_max)){  //CASE4: partial overlap (j=inner body), boundaries: inner=ri_min, outer=rj_max
        reb_integrator_hermes_autocalc_HSF_case_partial_overlap(mu, rj_min, rj_max, aj, ej, &vphi_max_j, &vr_max_j, ai, ei, ri_min, ri_max, &vphi_max_i, &vr_max_i);
    }
    //We calculate vrel_max the following way since it can be solved analytically. The correct way to find vrel_max is not easily
    //done (leads to a quartic soln). This estimate is guaranteed to be larger than the correct way, leading to a more conservative
    //estimate of min_dt_enc.
    double vrel_max2 = (vr_max_i+vr_max_j)*(vr_max_i+vr_max_j) + (vphi_max_i-vphi_max_j)*(vphi_max_i-vphi_max_j);
    if(vrel_max2 > 0.){
        double rhill_sum = ai*pow(r->particles[i].m/(3.*m0),1./3.) + aj*pow(r->particles[j].m/(3.*m0),1./3.);
        double dt_enc2 = rhill_sum*rhill_sum/vrel_max2;
        *min_dt_enc2 = MIN(*min_dt_enc2,dt_enc2);
    }
}

static void reb_integrator_hermes_autocalc_HSF_case_total_overlap(double mu, double ri_min, double ri_max, double ai, double ei, double* vphi_max_i, double* vr_max_i, double aj, double ej, double rj_min, double rj_max, double* vphi_max_j, double* vr_max_j){
    //body 'i' completely overlaps body 'j'
    const double sqrt_termj = sqrt(1.-ej*ej);
    const double sqrt_termi = sqrt(1.-ei*ei);
    const double ni = sqrt(mu/(ai*ai*ai));
    const double nj = sqrt(mu/(aj*aj*aj));
    *vphi_max_j = nj*aj*(1.+ej)/sqrt_termj;             //max phi velocity (j) is @ rj_min = aj*(1.-ej)
    *vphi_max_i = ni*ai*ai*sqrt_termi/(aj*(1.-ej));     //max phi velocity (i) is *also* @ rj_min = aj*(1.-ej), i.e. the smallest distance in the total overlap region
    *vr_max_j = nj*aj*ej/sqrt_termj;                    //max radial velocity (j) is @ rj = aj*(1.-ej^2.)
    *vr_max_i = ni*ai*ei/sqrt_termi;                    //max radial velocity (i) is @ ri = ai*(1.-ei^2.), *but check if this falls in overlap region*
    double loc_vr_max_i = ai*(1.-ei*ei);                //location of max radial velocity (i)
    if((loc_vr_max_i>rj_max)||(loc_vr_max_i<rj_min)){   //is loc_vr_max_i outside the overlap region?
        *vr_max_i *= sqrt(MAX( 1.-pow(loc_vr_max_i/(rj_min*ei)-1./ei,2.), 1.-pow(loc_vr_max_i/(rj_max*ei)-1./ei,2.) )); //if so, vr_max_i = MAX(vr(rj_min), vr(rj_max))
    }
}

static void reb_integrator_hermes_autocalc_HSF_case_partial_overlap(double mu, double ri_min, double ri_max, double ai, double ei, double* vphi_max_i, double* vr_max_i, double aj, double ej, double rj_min, double rj_max, double* vphi_max_j, double* vr_max_j){
    //body 'i' is interior to body 'j'
    const double sqrt_termj = sqrt(1.-ej*ej);
    const double sqrt_termi = sqrt(1.-ei*ei);
    const double ni = sqrt(mu/(ai*ai*ai));
    const double nj = sqrt(mu/(aj*aj*aj));
    *vphi_max_j = nj*aj*(1.+ej)/sqrt_termj;             //max phi velocity (j) is @ rj_min = aj*(1.-ej)
    *vphi_max_i = ni*ai*ai*sqrt_termi/(aj*(1.-ej));     //max phi velocity (i) is *also* @ rj_min = aj*(1.-ej), i.e. the smallest distance in the total overlap region
    *vr_max_j = nj*aj*ej/sqrt_termj;                    //max radial velocity (j) is @ rj = aj*(1.-ej^2.), *but check if this falls in overlap region*
    *vr_max_i = ni*ai*ei/sqrt_termi;                    //max radial velocity (i) is @ ri = ai*(1.-ei^2.), *but check if this falls in overlap region*
    double loc_vr_max_i = ai*(1.-ei*ei);                //location of max radial velocity (i) is @ ri = ai*(1.-ei^2.)
    double loc_vr_max_j = aj*(1.-ej*ej);                //location of max radial velocity (j) is @ rj = aj*(1.-ej^2.)
    if(loc_vr_max_i<rj_min){                            //is loc_vr_max_i outside the overlap region?
        *vr_max_i *= sqrt(1.-pow(loc_vr_max_i/(rj_min*ei)-1./ei,2.));  //if so, vr_max_i occurs at rj_min
    }
    if(loc_vr_max_j>ri_max){                            //is loc_vr_max_j outside the overlap region?
        *vr_max_j *= sqrt(1.-pow(loc_vr_max_j/(ri_max*ej)-1./ej,2.));  //if so, vr_max_i occurs at ri_max
    }
}

static void reb_integrator_hermes_apply_forces(const struct reb_simulation* r, double* a){
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
    struct reb_simulation* global = mini->ri_hermes.global;
    if (mini->testparticle_type){
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
    
    if(global->additional_forces){
        global->additional_forces(mini);
    }
}

static void reb_integrator_hermes_get_ae(struct reb_simulation* r, int index, double* a, double* e){
    struct reb_particle m0 = r->particles[0];
    struct reb_particle* particles = r->particles;
    const double G = r->G;
    const double mu = G*m0.m;
    const double muinv = 1./mu;
    
    struct reb_particle p = particles[index];
    const double dvx = p.vx-m0.vx;
    const double dvy = p.vy-m0.vy;
    const double dvz = p.vz-m0.vz;
    const double dx = p.x-m0.x;
    const double dy = p.y-m0.y;
    const double dz = p.z-m0.z;
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt(dx*dx + dy*dy + dz*dz);   //distance
    const double dinv = 1./d;
    const double vr = (dx*dvx + dy*dvy + dz*dvz)*dinv;
    const double ex = muinv*( (v2-mu*dinv)*dx - d*vr*dvx );
    const double ey = muinv*( (v2-mu*dinv)*dy - d*vr*dvy );
    const double ez = muinv*( (v2-mu*dinv)*dz - d*vr*dvz );
    
    *e = sqrt( ex*ex + ey*ey + ez*ez );    //eccentricity
    *a = -mu/(v2 - 2.*mu*dinv);            //semi-major axis
}
