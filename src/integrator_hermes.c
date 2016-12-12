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
static void reb_integrator_hermes_autocalc_HSF(struct reb_simulation* r);
static void reb_integrator_hermes_get_ae(struct reb_simulation* r, struct reb_particle com, int index, double* a, double* e);

void reb_integrator_hermes_part1(struct reb_simulation* r){
    r->gravity_ignore_terms = 0;
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

    // Determine HSF
    r->ri_hermes.current_hill_switch_factor = r->ri_hermes.hill_switch_factor;
    if(r->ri_hermes.adaptive_hill_switch_factor){
        reb_integrator_hermes_autocalc_HSF(r); // increases current_hill_switch_factor is needed
    }
    
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
        r->energy_offset += r->ri_hermes.mini->energy_offset;
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
            }
        }
    }
    if (global->ri_hermes.adaptive_hill_switch_factor==0 && global->ri_hermes.timestep_too_large_warning==0 && min_dt_enc2 < 16.*global->dt*global->dt){
        global->ri_hermes.timestep_too_large_warning = 1;
        reb_warning(global,"The timestep is likely too large. Close encounters might be missed. Decrease the timestep or increase the switching radius. This warning will appear only once.");
    }
}

//get min encounter time between overlapping orbits
static void reb_integrator_hermes_autocalc_HSF(struct reb_simulation* r){
    const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    struct reb_particle com = reb_get_com(r);
    const double mu = r->G*r->particles[0].m;
    struct reb_particle* particles = r->particles;
    double min_dt_enc2 = INFINITY;
    double m0 = particles[0].m;
    int* is_in_mini = r->ri_hermes.is_in_mini;
    for(int i=1;i<_N_active;i++){                                               //run over massive bodies
        double ep, ap;
        reb_integrator_hermes_get_ae(r, com, i, &ap, &ep);
        double rp_min = ap*(1.-ep);
        double rp_max = ap*(1.+ep);
        double np = sqrt(mu/(ap*ap*ap));
        for(int j=i+1;j<r->N;j++){                                              //run over massive + planetesimal bodies
            if(is_in_mini[j] == 1) continue;                                    //exclude bodies in mini from Auto HSF calc
            double e, a, n;
            reb_integrator_hermes_get_ae(r, com, j, &a, &e);
            double r_min = a*(1.-e);
            double r_max = a*(1.+e);
            double vphi_max_r=0., vr_max_r=0., global_max_r=0., sinf_max_r=0.;
            double vphi_max_rp=0., vr_max_rp=0., global_max_rp=0., sinf_max_rp=0.;
            if((rp_min<r_min)&&(rp_max>r_max)){         //CASE1: massive planet totally overlaps planetesimal
                n = sqrt(mu/(a*a*a));
                vphi_max_r = n*a*(1.+e)/sqrt(1.-e*e);                           //vphi_max is at r_min = a*(1.-e)
                vphi_max_rp = np*ap*ap*(1.-ep*ep)/(a*(1.-e)*sqrt(1.-ep*ep));    //vphi_max_rp @ r_min
                vr_max_r = n*a*e/sqrt(1.-e*e);                                  //vr_max_r is at r = a*(1.-e^2.)
                global_max_rp = ap*(1.-ep*ep);                                  //the distance corresponding to the global vr_max_rp
                if((global_max_rp>r_max)||(global_max_rp<r_min)){               //take max of boundaries (r_min and r_max)
                    sinf_max_rp = sqrt(MAX(1.-pow(global_max_rp/(r_min*ep)-1./ep,2.), 1.-pow(global_max_rp/(r_max*ep)-1./ep,2.)));
                    vr_max_rp = np*ap*ep/sqrt(1.-ep*ep) * sinf_max_rp;
                } else { vr_max_rp = np*ap*ep/sqrt(1.-ep*ep); }
            } else if((r_min<rp_min)&&(r_max>rp_max)){  //CASE2: planetesimal totally overlaps planet
                n = sqrt(mu/(a*a*a));
                vphi_max_rp = np*ap*(1.+ep)/sqrt(1.-ep*ep);
                vphi_max_r = n*a*a*(1.-e*e)/(ap*(1.-ep)*sqrt(1.-e*e));
                vr_max_rp = np*ap*ep/sqrt(1.-ep*ep);
                global_max_r = a*(1.-e*e);
                if((global_max_r>rp_max)||(global_max_r<rp_min)){               //take max of boundaries (rp_min and rp_max)
                    sinf_max_r = sqrt(MAX(1.-pow(global_max_r/(rp_min*e)-1./e,2.), 1.-pow(global_max_r/(rp_max*e)-1./e,2.)));
                    vr_max_r = n*a*e/sqrt(1.-e*e) * sinf_max_r;
                } else {vr_max_r = n*a*e/sqrt(1.-e*e);}
            } else if((rp_max>r_max)&&(r_max>rp_min)){  //CASE3: partial overlap (planetesimal=inner body), boundaries: inner=rp_min, outer=r_max
                n = sqrt(mu/(a*a*a));
                vphi_max_r = n*a*a*(1.-e*e)/(ap*(1.-ep)*sqrt(1.-e*e));
                vphi_max_rp = np*ap*(1.+ep)/sqrt(1.-ep*ep);
                global_max_r = a*(1.-e*e);
                if(global_max_r<rp_min){                                        //Since r_max is a minimum of vr, vr_max_r must be at rp_min
                    vr_max_r = n*a*e*sqrt((1.-pow(global_max_r/(rp_min*e)-1./e,2.))/(1.-e*e));
                } else {vr_max_r = n*a*e/sqrt(1.-e*e);}
                global_max_rp = ap*(1.-ep*ep);
                if(global_max_rp>r_max){                                        //Since rp_min is a minimum of vr, vr_max_rp must be at r_max
                    vr_max_rp = np*ap*ep*sqrt((1.-pow(global_max_rp/(r_max*ep)-1./ep,2.))/(1.-ep*ep));
                } else {vr_max_rp = np*ap*ep/sqrt(1.-ep*ep);}
            } else if((r_max>rp_max)&&(rp_max>r_min)){  //CASE4: partial overlap (planet=inner body), boundaries: inner=r_min, outer=rp_max
                n = sqrt(mu/(a*a*a));
                vphi_max_r = n*a*(1.+e)/sqrt(1.-e*e);
                vphi_max_rp = np*ap*ap*(1.-ep*ep)/(a*(1.-e)*sqrt(1.-ep*ep));
                global_max_r = a*(1.-e*e);
                if(global_max_r>rp_max){                                        //Since r_min is a minimum of vr, vr_max_r must be at rp_max
                    vr_max_r = n*a*e*sqrt((1.-pow(global_max_r/(rp_max*e)-1./e,2.))/(1.-e*e));
                } else {vr_max_r = n*a*e/sqrt(1.-e*e);}
                global_max_rp = ap*(1.-ep*ep);
                if(global_max_rp<r_min){                                        //Since rp_max is a minimum of vr, vr_max_rp must be at r_min
                    vr_max_rp = np*ap*ep*sqrt((1.-pow(global_max_rp/(r_min*ep)-1./ep,2.))/(1.-ep*ep));
                } else {vr_max_rp = np*ap*ep/sqrt(1.-ep*ep);}
            }
            //We calculate vrel_max the following way since it can be solved analytically. The correct way to find vrel_max is not easily
            //done (leads to a quartic soln). This estimate is guaranteed to be larger than the correct way, leading to a more conservative
            //estimate of min_dt_enc.
            double vrel_max2 = (vr_max_rp+vr_max_r)*(vr_max_rp+vr_max_r) + (vphi_max_rp-vphi_max_r)*(vphi_max_rp-vphi_max_r);
            if(vrel_max2 > 0.){
                double rhill_sum = ap*pow(particles[i].m/(3.*m0),1./3.) + a*pow(particles[j].m/(3.*m0),1./3.);
                double dt_enc2 = rhill_sum*rhill_sum/vrel_max2;
                min_dt_enc2 = MIN(min_dt_enc2,dt_enc2);
            }
        }
    }
    
    if(min_dt_enc2 < INFINITY){
        double dt2 = 16.*r->dt*r->dt;                                           //Factor of 4:hill sphere > 4*length scales for wiggle room
        double HSF_new = sqrt(dt2/min_dt_enc2);
        double base = 1.25;
        double exp = ceilf(log10(HSF_new)/log10(base));                         //round HSF up to nearest multiple of 1.25
        HSF_new = pow(base,exp);
        r->ri_hermes.current_hill_switch_factor = MAX(r->ri_hermes.current_hill_switch_factor, HSF_new); // Increase HSF if needed
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

static void reb_integrator_hermes_get_ae(struct reb_simulation* r, struct reb_particle com, int index, double* a, double* e){
    const double G = r->G;
    const double mu = G*r->particles[0].m;
    const double muinv = 1./mu;
    struct reb_particle* particles = r->particles;
    
    struct reb_particle p = particles[index];
    const double dvx = p.vx-com.vx;
    const double dvy = p.vy-com.vy;
    const double dvz = p.vz-com.vz;
    const double dx = p.x-com.x;
    const double dy = p.y-com.y;
    const double dz = p.z-com.z;
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
