/**
 * @file    integrator_mercurius.c
 * @brief   MERCURIUS, a modified version of John Chambers' MERCURY algorithm
 *          using the IAS15 integrator and WHFast. It works with planet-planry
 *          collisions, test particles, and additional forces.
 * @author  Hanno Rein, Dan Tamayo
 * 
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein, Dan Tamayo
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
#include "integrator.h"
#include "gravity.h"
#include "integrator_mercurius.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function used by the Mercury integrator.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return 10.*(y*y*y) - 15.*(y*y*y*y) + 6.*(y*y*y*y*y);
    }
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}

double reb_integrator_mercurius_L_infinity(const struct reb_simulation* const r, double d, double dcrit){
    // Infinitely differentiable function.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
}


void reb_integrator_mercurius_inertial_to_dh(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_vec3d com_pos = {0};
    struct reb_vec3d com_vel = {0};
    double mtot = 0.;
    const int N_active = r->N_active==-1?r->N:r->N_active;
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
    r->ri_mercurius.com_pos = com_pos;
    r->ri_mercurius.com_vel = com_vel;
}

void reb_integrator_mercurius_dh_to_inertial(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle temp = {0};
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
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
    particles[0].x = r->ri_mercurius.com_pos.x - temp.x; 
    particles[0].y = r->ri_mercurius.com_pos.y - temp.y; 
    particles[0].z = r->ri_mercurius.com_pos.z - temp.z; 

    for (int i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += r->ri_mercurius.com_vel.x;
        particles[i].vy += r->ri_mercurius.com_vel.y;
        particles[i].vz += r->ri_mercurius.com_vel.z;
    }
    particles[0].vx = r->ri_mercurius.com_vel.x - temp.vx; 
    particles[0].vy = r->ri_mercurius.com_vel.y - temp.vy; 
    particles[0].vz = r->ri_mercurius.com_vel.z - temp.vz; 
}


static void reb_mercurius_encounter_predict(struct reb_simulation* const r){
    // This function predicts close encounters during the timestep
    // It makes use of the old and new position and velocities obtained
    // after the Kepler step.
    struct reb_simulation_integrator_mercurius* rim = &(r->ri_mercurius);
    struct reb_particle* const particles = r->particles;
    struct reb_particle* const particles_backup = rim->particles_backup;
    const double* const dcrit = rim->dcrit;
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const double dt = r->dt;
    rim->encounterN = 1;
    rim->encounter_map[0] = 1;
    for (int i=1; i<N; i++){
        rim->encounter_map[i] = 0;
    }
    for (int i=0; i<N_active; i++){
        for (int j=i+1; j<N; j++){
            const double dxn = particles[i].x - particles[j].x;
            const double dyn = particles[i].y - particles[j].y;
            const double dzn = particles[i].z - particles[j].z;
            const double dvxn = particles[i].vx - particles[j].vx;
            const double dvyn = particles[i].vy - particles[j].vy;
            const double dvzn = particles[i].vz - particles[j].vz;
            const double rn = (dxn*dxn + dyn*dyn + dzn*dzn);
            const double dxo = particles_backup[i].x - particles_backup[j].x;
            const double dyo = particles_backup[i].y - particles_backup[j].y;
            const double dzo = particles_backup[i].z - particles_backup[j].z;
            const double dvxo = particles_backup[i].vx - particles_backup[j].vx;
            const double dvyo = particles_backup[i].vy - particles_backup[j].vy;
            const double dvzo = particles_backup[i].vz - particles_backup[j].vz;
            const double ro = (dxo*dxo + dyo*dyo + dzo*dzo);

            const double drndt = (dxn*dvxn+dyn*dvyn+dzn*dvzn)*2.;
            const double drodt = (dxo*dvxo+dyo*dvyo+dzo*dvzo)*2.;

            const double a = 6.*(ro-rn)+3.*dt*(drodt+drndt); 
            const double b = 6.*(rn-ro)-2.*dt*(2.*drodt+drndt); 
            const double c = dt*drodt; 

            double rmin = MIN(rn,ro);

            const double s = b*b-4.*a*c;
            const double sr = sqrt(MAX(0.,s));
            const double tmin1 = (-b + sr)/(2.*a); 
            const double tmin2 = (-b - sr)/(2.*a); 
            if (tmin1>0. && tmin1<1.){
                const double rmin1 = (1.-tmin1)*(1.-tmin1)*(1.+2.*tmin1)*ro
                                     + tmin1*tmin1*(3.-2.*tmin1)*rn
                                     + tmin1*(1.-tmin1)*(1.-tmin1)*dt*drodt
                                     - tmin1*tmin1*(1.-tmin1)*dt*drndt;
                rmin = MIN(MAX(rmin1,0.),rmin);
            }
            if (tmin2>0. && tmin2<1.){
                const double rmin2 = (1.-tmin2)*(1.-tmin2)*(1.+2.*tmin2)*ro
                                     + tmin2*tmin2*(3.-2.*tmin2)*rn
                                     + tmin2*(1.-tmin2)*(1.-tmin2)*dt*drodt
                                     - tmin2*tmin2*(1.-tmin2)*dt*drndt;
                rmin = MIN(MAX(rmin2,0.),rmin);
            }

            if (sqrt(rmin)< 1.1*MAX(dcrit[i],dcrit[j])){
                if (rim->encounter_map[i]==0){
                    rim->encounter_map[i] = i;
                    rim->encounterN++;
                }
                if (rim->encounter_map[j]==0){
                    rim->encounter_map[j] = j;
                    rim->encounterN++;
                }
            }
        }
    }
}
    
void reb_integrator_mercurius_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_mercurius_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    double px=0., py=0., pz=0.;
    for (int i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m; 
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px /= r->particles[0].m;
    py /= r->particles[0].m;
    pz /= r->particles[0].m;
    for (int i=1;i<N;i++){
        particles[i].x += dt*px;
        particles[i].y += dt*py;
        particles[i].z += dt*pz;
    }
}

void reb_integrator_mercurius_com_step(struct reb_simulation* const r, double dt){
    r->ri_mercurius.com_pos.x += dt*r->ri_mercurius.com_vel.x;
    r->ri_mercurius.com_pos.y += dt*r->ri_mercurius.com_vel.y;
    r->ri_mercurius.com_pos.z += dt*r->ri_mercurius.com_vel.z;
}

void reb_integrator_mercurius_kepler_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        reb_whfast_kepler_solver(r,particles,r->G*particles[0].m,i,dt); // in dh
    }
}

static void reb_mercurius_encounter_step(struct reb_simulation* const r, const double _dt){
    // Only particles having a close encounter are integrated by IAS15.
    struct reb_simulation_integrator_mercurius* rim = &(r->ri_mercurius);
    if (rim->encounterN<2){
        return; // If there are no particles (other than the star) having a close encounter, then there is nothing to do.
    }

    int i_enc = 0;
    rim->encounterNactive = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(rim->encounter_map[i]){  
            r->particles[i] = rim->particles_backup[i]; // use coordinates before whfast step
            rim->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
                rim->encounterNactive++;
            }
        }
    }

    rim->mode = 1;
    
    // run
    const double old_dt = r->dt;
    const double old_t = r->t;
    double t_needed = r->t + _dt; 
        
    reb_integrator_ias15_reset(r);
    
    r->dt = 0.0001*_dt; // start with a small timestep.
    
    while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 ){
        struct reb_particle star = r->particles[0]; // backup velocity
        r->particles[0].vx = 0; // star does not move in dh 
        r->particles[0].vy = 0;
        r->particles[0].vz = 0;
        reb_update_acceleration(r);
        reb_integrator_ias15_part2(r);
        r->particles[0].vx = star.vx; // restore every timestep for collisions
        r->particles[0].vy = star.vy;
        r->particles[0].vz = star.vz;
        
        if (r->t+r->dt >  t_needed){
            r->dt = t_needed-r->t;
        }

        // Search and resolve collisions
        reb_collision_search(r);

        // Do any additional post_timestep_modifications.
        // Note: post_timestep_modifications is called here but also
        // at the end of the full timestep. The function thus needs
        // to be implemented with care as not to do the same 
        // modification multiple times. To do that, check the value of
        // r->ri_mercurius.mode
        if (r->post_timestep_modifications){
            r->post_timestep_modifications(r);
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

    // Reset constant for global particles
    r->t = old_t;
    r->dt = old_dt;
    rim->mode = 0;

}

void reb_integrator_mercurius_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }
    
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    const int N = r->N;
    
    if (rim->dcrit_allocatedN<N){
        // Need to safe these arrays in SimulationArchive
        rim->dcrit              = realloc(rim->dcrit, sizeof(double)*N);
        rim->dcrit_allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_dcrit_this_timestep        = 1;
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
        rim->recalculate_coordinates_this_timestep = 1;
    }
    if (rim->allocatedN<N){
        // These arrays are only used within one timestep. 
        // Can be recreated without loosing bit-wise reproducibility
        rim->particles_backup   = realloc(rim->particles_backup,sizeof(struct reb_particle)*N);
        rim->encounter_map      = realloc(rim->encounter_map,sizeof(int)*N);
        rim->allocatedN = N;
    }
    if (rim->safe_mode || rim->recalculate_coordinates_this_timestep){
        if (rim->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r);
            reb_warning(r,"MERCURIUS: Recalculating heliocentric coordinates but coordinates were not synchronized before.");
        }
        reb_integrator_mercurius_inertial_to_dh(r);
        rim->recalculate_coordinates_this_timestep = 0;
    }

    if (rim->recalculate_dcrit_this_timestep){
        rim->recalculate_dcrit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r);
            reb_warning(r,"MERCURIUS: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        rim->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        const double m0 = r->particles[0].m;
        for (int i=1;i<N;i++){
            const double dx  = r->particles[i].x;  // in dh
            const double dy  = r->particles[i].y;
            const double dz  = r->particles[i].z;
            const double dvx = r->particles[i].vx - r->particles[0].vx; 
            const double dvy = r->particles[i].vy - r->particles[0].vy; 
            const double dvz = r->particles[i].vz - r->particles[0].vz; 
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);
            const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;

            const double GM = r->G*(m0+r->particles[i].m);
            const double a = GM*_r / (2.*GM - _r*v2);
            const double vc = sqrt(GM/fabs(a));
            double dcrit = 0;
            // Criteria 1: average velocity
            dcrit = MAX(dcrit, vc*0.4*r->dt);
            // Criteria 2: current velocity
            dcrit = MAX(dcrit, sqrt(v2)*0.4*r->dt);
            // Criteria 3: Hill radius
            dcrit = MAX(dcrit, rim->hillfac*a*pow(r->particles[i].m/(3.*r->particles[0].m),1./3.));
            // Criteria 4: physical radius
            dcrit = MAX(dcrit, 2.*r->particles[i].r);

            rim->dcrit[i] = dcrit;
        }
    }
    
    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }
    
    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_MERCURIUS){
        reb_warning(r,"Mercurius has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_MERCURIUS;
    rim->mode = 0;
    
    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurius_L_mercury;
    }
}

void reb_integrator_mercurius_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    const int N = r->N;
   
    if (rim->is_synchronized){
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);
    }else{
        reb_integrator_mercurius_interaction_step(r,r->dt);
    }
    reb_integrator_mercurius_jump_step(r,r->dt/2.);
    reb_integrator_mercurius_com_step(r,r->dt); 
    
    // Make copy of particles before the kepler step.
    // Then evolve all particles in kepler step.
    // Result will be used in encounter prediction.
    // Particles having a close encounter will be overwritten 
    // later by encounter step.
    memcpy(rim->particles_backup,r->particles,N*sizeof(struct reb_particle)); 
    reb_integrator_mercurius_kepler_step(r,r->dt);
    
    reb_mercurius_encounter_predict(r);
   
    reb_mercurius_encounter_step(r,r->dt);
    
    reb_integrator_mercurius_jump_step(r,r->dt/2.);
        
    rim->is_synchronized = 0;
    if (rim->safe_mode){
        reb_integrator_mercurius_synchronize(r);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

void reb_integrator_mercurius_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_mercurius* const rim = &(r->ri_mercurius);
    if (rim->is_synchronized == 0){
        r->gravity = REB_GRAVITY_MERCURIUS; // needed here again for SimulationArchive
        rim->mode = 0;
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurius_L_mercury;
        }
        reb_update_acceleration(r);
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);
        
        reb_integrator_mercurius_dh_to_inertial(r);

        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurius_reset(struct reb_simulation* r){
    r->ri_mercurius.L = NULL;
    r->ri_mercurius.mode = 0;
    r->ri_mercurius.encounterN = 0;
    r->ri_mercurius.encounterNactive = 0;
    r->ri_mercurius.hillfac = 3;
    r->ri_mercurius.recalculate_coordinates_this_timestep = 0;
    // Internal arrays (only used within one timestep)
    free(r->ri_mercurius.particles_backup);
    r->ri_mercurius.particles_backup = NULL;
    free(r->ri_mercurius.particles_backup_additionalforces);
    r->ri_mercurius.particles_backup_additionalforces = NULL;
    free(r->ri_mercurius.encounter_map);
    r->ri_mercurius.encounter_map = NULL;
    r->ri_mercurius.allocatedN = 0;
    r->ri_mercurius.allocatedN_additionalforces = 0;
    // dcrit array
    free(r->ri_mercurius.dcrit);
    r->ri_mercurius.dcrit = NULL;
    r->ri_mercurius.dcrit_allocatedN = 0;
}

