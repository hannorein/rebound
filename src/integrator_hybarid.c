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
void reb_integrator_hybarid_additional_forces_mini(struct reb_simulation* mini);
void mini_check_for_collision(struct reb_simulation* mini);

double E0;

void reb_integrator_hybarid_part1(struct reb_simulation* r){
	const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    if (r->ri_hybarid.mini == NULL){
        r->ri_hybarid.mini = reb_create_simulation();
        r->ri_hybarid.mini->usleep = -1; // Disable visualiation
        r->ri_hybarid.mini->integrator = REB_INTEGRATOR_IAS15;
        r->ri_hybarid.mini->additional_forces = reb_integrator_hybarid_additional_forces_mini;
        r->ri_hybarid.mini->ri_hybarid.global = r;
        r->ri_hybarid.mini->testparticle_type = r->testparticle_type;
        r->ri_hybarid.mini->heartbeat = mini_check_for_collision;
        //r->ri_hybarid.mini->softening = 1e-6;
        //r->ri_hybarid.mini->ri_ias15.epsilon = 1e-8;  //speeds up ias and hybarid immensely
        E0 = reb_tools_energy(r);
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

    //keep this after check_for_encounter - then if particle is removed, no need to edit particles_prev
    if (r->testparticle_type){
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
	// Do nothing.
    reb_integrator_whfast_reset(r);
}

static void reb_integrator_hybarid_check_for_encounter(struct reb_simulation* r){
    struct reb_simulation* mini = r->ri_hybarid.mini;
	const int _N_active = ((r->N_active==-1)?r->N:r->N_active) - r->N_var;
    struct reb_particle* global = r->particles;
    struct reb_particle p0 = global[0];
    double ejectiondistance2 = 100;     //temporary hardcoded value.
    double HSR = r->ri_hybarid.switch_ratio;
    double minr=100; double max_vr = 1e-6;  int piid = 0; int pjid=0;//AS temp
    for (int i=0; i<_N_active; i++){
        struct reb_particle* pi = &(global[i]);
        double rhi;
        if(i==0) rhi = 100*p0.r*p0.r; else{
            const double dxi = p0.x - pi->x;
            const double dyi = p0.y - pi->y;
            const double dzi = p0.z - pi->z;
            const double r0i2 = dxi*dxi + dyi*dyi + dzi*dzi;
            rhi = r0i2*pow(pi->m/(p0.m*3.),2./3.);
        }
        for(int j=i+1;j<r->N;j++){
            struct reb_particle pj = global[j];
            
            const double dxj = p0.x - pj.x;
            const double dyj = p0.y - pj.y;
            const double dzj = p0.z - pj.z;
            const double r0j2 = dxj*dxj + dyj*dyj + dzj*dzj;
            const double rhj = r0j2*pow(pj.m/(p0.m*3.),2./3.);  //this is calculated for each massive planet but only needs to be once.
            
            const double dx = pi->x - pj.x;
            const double dy = pi->y - pj.y;
            const double dz = pi->z - pj.z;
            const double rij2 = dx*dx + dy*dy + dz*dz;
            const double ratio = rij2/(rhi+rhj);    //(p-p distance/Hill radii)^2

            if(ratio < HSR){
                r->ri_hybarid.mini_active = 1;
                if (j>=_N_active && r->ri_hybarid.is_in_mini[j] ==0){//make sure not already added
                    reb_add(mini,pj);
                    r->ri_hybarid.is_in_mini[j] = 1;
                    if (r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax){
                        while(r->ri_hybarid.encounter_index_N>=r->ri_hybarid.encounter_index_Nmax) r->ri_hybarid.encounter_index_Nmax += 32;
                        r->ri_hybarid.encounter_index = realloc(r->ri_hybarid.encounter_index,r->ri_hybarid.encounter_index_Nmax*sizeof(int));
                    }
                    r->ri_hybarid.encounter_index[r->ri_hybarid.encounter_index_N] = j;
                    r->ri_hybarid.encounter_index_N++;
                }
            } else if (i==0 && r0j2 > ejectiondistance2){
                double Ei = reb_tools_energy(r);
                reb_remove(r,j,1);
                double Ef = reb_tools_energy(r);
                r->ri_hybarid.dE_offset += Ei - Ef;
                printf("\n\tParticle %d ejected from system at t=%f, E=%e\n",pj.id,r->t,fabs((Ef+r->ri_hybarid.dE_offset-E0)/E0));
                j--;    //re-try iteration j since j+1 is now j but hasn't been checked.
            }
            if(r->t > 15638.38 && r->t < 15638.49){
                double vx = pi->vx - pj.vx;
                double vy = pi->vy - pj.vy;
                double vz = pi->vz - pj.vz;
                double vrel = sqrt(vx*vx + vy*vy + vz*vz);
                double rr = sqrt(rij2);
                double val = r->dt*vrel/rr;
                if(rr < minr) minr = rr;
                if(val > max_vr){
                    max_vr = val;
                    piid = pi->id;
                    pjid = pj.id;
                }
                double E = reb_tools_energy(r) + r->ri_hybarid.dE_offset;
                double dE = fabs((E-E0)/E0);
                
                FILE *append;
                append = fopen("output/debug.txt", "a");
                fprintf(append, "%.16f,%.16f,%d,%d,%f,%f,%d,%d,%d,%d,%d,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%f\n",r->t,dE,i,j,minr,max_vr,piid,pjid,r->N,mini->N,r->ri_hybarid.encounter_index_N,r->ri_hybarid.mini_active,rhi,r0j2,rhj,rij2,ratio,HSR);
                fclose(append);
            }
        }
    }
    //if(r->t > 15500 && r->t < 16005){
    //}
}

void reb_integrator_hybarid_additional_forces_mini(struct reb_simulation* mini){
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

//check for collisions in mini each heartbeat
void mini_check_for_collision(struct reb_simulation* mini){
    struct reb_simulation* r = mini->ri_hybarid.global;
    struct reb_particle* particles = mini->particles;
    int N_active = mini->N_active;
    double dtmini = mini->dt;
    for(int i=0;i<N_active;i++){
        struct reb_particle* pi = &(particles[i]);
        for(int j=N_active;j<mini->N;j++){
            struct reb_particle* pj = &(particles[j]);
            double radius2 = (pi->r+pj->r)*(pi->r+pj->r);
            double dvx = pi->vx - pj->vx;
            double dvy = pi->vy - pj->vy;
            double dvz = pi->vz - pj->vz;
            double dx = pj->x - pi->x;
            double dy = pj->y - pi->y;
            double dz = pj->z - pi->z;
            
            double tmin = (dx*dvx + dy*dvy + dz*dvz)/(dvx*dvx + dvy*dvy + dvz*dvz);
            double dmin2 = radius2*2.;
            if(tmin < dtmini) dmin2 = (dx - dvx*tmin)*(dx - dvx*tmin) + (dy - dvy*tmin)*(dy - dvy*tmin) + (dz - dvz*tmin)*(dz - dvz*tmin);
            
            double dstart2 = dx*dx + dy*dy + dz*dz;
            double dend2 = (dx - dvx*dtmini)*(dx - dvx*dtmini) + (dy - dvy*dtmini)*(dy - dvy*dtmini) + (dz - dvz*dtmini)*(dz - dvz*dtmini);
            if(dmin2 <= radius2 || dstart2 <= radius2 || dend2 <= radius2){
                //pj->lastcollision = mini->t;
                
                double invmass = 1.0/(pi->m + pj->m);
                double Ei = reb_tools_energy(mini);
                
                pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
                pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
                pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
                pi->m += pj->m;
                
                reb_remove(mini,j,1);  //remove mini
                
                double Ef = reb_tools_energy(mini);
                r->ri_hybarid.dE_offset += Ei - Ef;
                printf("\n\tParticle %d collided with body %d from system at t=%f.\n",i,r->ri_hybarid.encounter_index[j],r->t);
                
                //remove from global and update global arrays
                int globalj = r->ri_hybarid.encounter_index[j];
                reb_remove(r,globalj,1);
                
                for(int k=globalj;k<r->N;k++){
                    r->ri_hybarid.particles_prev[k] = r->ri_hybarid.particles_prev[k+1];
                    r->ri_hybarid.is_in_mini[k] = r->ri_hybarid.is_in_mini[k+1];
                }
                r->ri_hybarid.encounter_index_N--;
                for(int k=j;k<r->ri_hybarid.encounter_index_N;k++) r->ri_hybarid.encounter_index[k] = r->ri_hybarid.encounter_index[k+1];
                for(int k=N_active;k<r->ri_hybarid.encounter_index_N;k++){
                    if(r->ri_hybarid.encounter_index[k] > globalj) r->ri_hybarid.encounter_index[k]--; //1 fewer particles in index now
                }
                
                j--;    //re-try iteration j since j+1 is now j but hasn't been checked.
            }
        }
    }
}
