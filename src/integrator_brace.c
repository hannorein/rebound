/**
 * @file    integrator_brace.c
 * @brief   BRACE
 * @author  Tiger Lu, Hanno Rein
 *
 * @section LICENSE
 * Copyright (c) 2023 Tiger Lu
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
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_brace.h"
#include "integrator_whfast.h"
#include "integrator_bs.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

int reb_integrator_brace_switch_default(struct reb_simulation* const r, const unsigned int i, const unsigned int j){
    // Returns 1 for close encounter between i and j, 0 otherwise
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    const double h2 = r->dt/2.;
    
    const double dxi  = r->particles[i].x;
    const double dyi  = r->particles[i].y;
    const double dzi  = r->particles[i].z;

    const double dxj  = r->particles[j].x;
    const double dyj  = r->particles[j].y;
    const double dzj  = r->particles[j].z;

    const double dx = dxi - dxj;
    const double dy = dyi - dyj;
    const double dz = dzi - dzj;
    const double rp = dx*dx + dy*dy + dz*dz;
    
    double dcriti6 = 0.0;
    double dcritj6 = 0.0;

    const double m0 = r->particles[0].m;
    
    // Check central body for physical radius ONLY
    if (i == 0 && r->particles[i].r != 0){
	const double rs = r->particles[0].r;
	dcriti6 = rs*rs*rs*rs*rs*rs;
    }

    else if (r->particles[i].m != 0){
        const double di2 = dxi*dxi + dyi*dyi + dzi*dzi;
        const double mr = r->particles[i].m/(3.*m0);
        dcriti6 = di2*di2*di2*mr*mr;
    }

    if (r->particles[j].m != 0){
        const double dj2 = dxj*dxj + dyj*dyj + dzj*dzj;
        const double mr = r->particles[j].m/(3.*m0);
        dcritj6 = dj2*dj2*dj2*mr*mr;
    }

    double r_crit_hill2 = ri_brace->r_crit_hill*ri_brace->r_crit_hill;
    double dcritmax6 = r_crit_hill2 * r_crit_hill2 * r_crit_hill2 * MAX(dcriti6,dcritj6);

    if (rp*rp*rp < dcritmax6) return 1;
    
    const double dvx  = r->particles[i].vx - r->particles[j].vx;
    const double dvy  = r->particles[i].vy - r->particles[j].vy;
    const double dvz  = r->particles[i].vz - r->particles[j].vz;
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    
    const double qv = dx*dvx + dy*dvy + dz*dvz;
    int d;

    if (qv == 0.0){ // Small
        // minimum is at present, which is already checked for
	return 0;
    }
    else if (qv < 0){
        d = 1; 
    }
    else{
        d = -1;
    }

    double dmin2;
    double tmin = -d*qv/v2;
    if (tmin < h2){
	// minimum is in the window
	dmin2 = rp - qv*qv/v2;
    }
    else{
	dmin2 = rp + 2*d*qv*h2 + v2*h2*h2;
    }

    return dmin2*dmin2*dmin2 < dcritmax6;
}

int reb_integrator_brace_switch_peri_default(struct reb_simulation* const r, const unsigned int j){
    // Following Pham et al (2024)
    const struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    double GM = r->G*r->particles[0].m; // Not sure if this is the right mass to use.
    
    double x = r->particles[j].x;
    double y = r->particles[j].y;
    double z = r->particles[j].z;
    double d2 = x*x + y*y + z*z;
    double d = sqrt(d2);

    // first derivative
    double dx = r->particles[j].vx;
    double dy = r->particles[j].vy;
    double dz = r->particles[j].vz;

    // second derivative
    double prefact2 = -GM/(d2*d);
    double ddx = prefact2*x;
    double ddy = prefact2*y;
    double ddz = prefact2*z;
    // need sqrt for this one...
    double dd = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);

    // third derivative
    double prefact3 = GM/(d2*d2*d);
    double dddx = prefact3*(-dx*(y*y+z*z) + 2.*x*x*dx+3.*x*(y*dy+z*dz));
    double dddy = prefact3*(-dy*(x*x+z*z) + 2.*y*y*dy+3.*y*(x*dx+z*dz));
    double dddz = prefact3*(-dz*(x*x+y*y) + 2.*z*z*dz+3.*z*(x*dx+y*dy));

    double ddd2 = dddx*dddx + dddy*dddy + dddz*dddz;

    // fourth derivative
    double prefact4 = GM/(d2*d2*d2*d);
    double ddddx = prefact4* (d2 * (-ddx*(y*y+z*z) + 2.*x*x*ddx + dx*(y*dy + z*dz) + x*(4.*dx*dx + 3.*(y*ddy + dy*dy + z*ddz + dz*dz ))) - 5.*(x*dx+y*dy+z*dz)*(-dx*(y*y+z*z)+2.*x*x*dx + 3.*x*(y*dy+z*dz)));
    double ddddy = prefact4* (d2 * (-ddy*(x*x+z*z) + 2.*y*y*ddy + dy*(x*dx + z*dz) + y*(4.*dy*dy + 3.*(x*ddx + dx*dx + z*ddz + dz*dz ))) - 5.*(y*dy+x*dx+z*dz)*(-dy*(x*x+z*z)+2.*y*y*dy + 3.*y*(x*dx+z*dz)));
    double ddddz = prefact4* (d2 * (-ddz*(y*y+x*x) + 2.*z*z*ddz + dz*(y*dy + x*dx) + z*(4.*dz*dz + 3.*(y*ddy + dy*dy + x*ddx + dx*dx ))) - 5.*(z*dz+y*dy+x*dx)*(-dz*(y*y+x*x)+2.*z*z*dz + 3.*z*(y*dy+x*dx)));
    double dddd = sqrt(ddddx*ddddx + ddddy*ddddy + ddddz*ddddz);

    double tau_prs2 = 2.*dd*dd/(ddd2+dd*dddd); // Eq 16
    double dt_prs2 = ri_brace->peri_crit_eta * ri_brace->peri_crit_eta * tau_prs2;

    if (r->dt * r->dt > dt_prs2){
        return 1;
    }else{
        return 0;
    }
}

int reb_integrator_brace_switch_peri_none(struct reb_simulation* const r, const unsigned int j){
    // No pericenter flags
    return 0;
}

void reb_integrator_brace_inertial_to_barycentric(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_vec3d com_pos = {0};
    struct reb_vec3d com_vel = {0};
    double mtot = 0.;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;
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
    for (int i=0;i<N;i++){
        particles[i].x -= com_pos.x;
        particles[i].y -= com_pos.y;
        particles[i].z -= com_pos.z;
        particles[i].vx -= com_vel.x;
        particles[i].vy -= com_vel.y;
        particles[i].vz -= com_vel.z;
    }
    r->ri_brace.com_pos = com_pos;
    r->ri_brace.com_vel = com_vel;
}

void reb_integrator_brace_barycentric_to_inertial(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    
    for (int i=0;i<N;i++){
        particles[i].x += r->ri_brace.com_pos.x;
        particles[i].y += r->ri_brace.com_pos.y;
        particles[i].z += r->ri_brace.com_pos.z;
        particles[i].vx += r->ri_brace.com_vel.x;
        particles[i].vy += r->ri_brace.com_vel.y;
        particles[i].vz += r->ri_brace.com_vel.z;
    }
}

void reb_integrator_brace_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    r->ri_brace.mode = REB_BRACE_MODE_INTERACTION;
    reb_simulation_update_acceleration(r);
    for (int i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_brace_com_step(struct reb_simulation* const r, double dt){
    r->ri_brace.com_pos.x += dt*r->ri_brace.com_vel.x;
    r->ri_brace.com_pos.y += dt*r->ri_brace.com_vel.y;
    r->ri_brace.com_pos.z += dt*r->ri_brace.com_vel.z;
}

void reb_integrator_brace_whfast_step(struct reb_simulation* const r, double dt){
    //struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    double m = r->particles[0].m;
    for (int i=1;i<N;i++){
        m += r->particles[i].m;
    }
    for (int i=1;i<N;i++){
        reb_whfast_kepler_solver(r,r->particles,r->G*m,i,dt);
    }
}

void reb_integrator_brace_update_particles(struct reb_simulation* r, const double* y){
    int N = r->ri_brace.encounter_N;
    int* map = r->ri_brace.encounter_map;

    for (int i=0; i<N; i++){
        int mi = map[i];
        struct reb_particle* const p = &(r->particles[mi]);
        p->x  = y[i*6+0];
        p->y  = y[i*6+1];
        p->z  = y[i*6+2];
        p->vx = y[i*6+3];
        p->vy = y[i*6+4];
        p->vz = y[i*6+5];
    }
}

void reb_integrator_brace_nbody_derivatives_barycentric(struct reb_ode* ode, double* const yDot, const double* const y, double const t){
    struct reb_simulation* const r = ode->r;
    // BRACE always needs this to ensure the right Hamiltonian is evolved
    reb_integrator_brace_update_particles(r, y);
    reb_simulation_update_acceleration(r);

    int* map = r->ri_brace.encounter_map;
    int N = r->ri_brace.encounter_N;

    if (map==NULL){
        reb_simulation_error(r, "Cannot access BRACE map from BS.");
        return;
    }

    for (int i=0; i<N; i++){
        int mi = map[i];
        const struct reb_particle p = r->particles[mi];
        yDot[i*6+0] = p.vx; // Already checked for current_L
        yDot[i*6+1] = p.vy;
        yDot[i*6+2] = p.vz;
        yDot[i*6+3] = p.ax;
        yDot[i*6+4] = p.ay;
        yDot[i*6+5] = p.az;
    }
}


void reb_integrator_brace_bs_step(struct reb_simulation* const r, double dt){
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);

    if (ri_brace->encounter_N < 1){
        // No close encounters, skip
        return;
    }
    //printf("encounter: %d\n", ri_brace->encounter_N);

    int i_enc = 0;
    ri_brace->encounter_N_active = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(ri_brace->encounter_map[i]){
            r->particles[i] = ri_brace->particles_backup_kepler[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            ri_brace->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
                ri_brace->encounter_N_active++;
            }
        }
    }

    ri_brace->mode = REB_BRACE_MODE_KEPLER;
    
    // Only Partial BS uses this step 
    // run
    const double old_dt = r->dt;
    const double old_t = r->t;
    const double t_needed = r->t + dt;
    reb_integrator_bs_reset(r);

    // Temporarily remove all odes for BS step
    struct reb_ode** odes_backup = r->odes;
    int N_allocated_odes_backup = r->N_allocated_odes;
    int N_odes_backup = r->N_odes;
    r->odes = NULL;
    r->N_allocated_odes = 0;
    r->N_odes = 0;

    // Temporarily add new nbody ode for BS step
    struct reb_ode* nbody_ode = reb_ode_create(r, ri_brace->encounter_N*3*2);
    nbody_ode->derivatives = reb_integrator_brace_nbody_derivatives_barycentric;
    nbody_ode->needs_nbody = 0;

    // TODO: Support backwards integrations
    while(r->t < t_needed && fabs(dt/old_dt)>1e-14 ){
        double* y = nbody_ode->y;

        // In case of overshoot
        if (r->t + dt >  t_needed){
            dt = t_needed - r->t;
        }

        for (unsigned int i=0; i<ri_brace->encounter_N; i++){
            const int mi = ri_brace->encounter_map[i];
            const struct reb_particle p = r->particles[mi];
            y[i*6+0] = p.x;
            y[i*6+1] = p.y;
            y[i*6+2] = p.z;
            y[i*6+3] = p.vx;
            y[i*6+4] = p.vy;
            y[i*6+5] = p.vz;
        }

        int success = reb_integrator_bs_step(r, dt);
        if (success){
            r->t += dt;
        }
        dt = r->ri_bs.dt_proposed;
        reb_integrator_brace_update_particles(r, nbody_ode->y);

        reb_collision_search(r);
        if (r->N_allocated_collisions) ri_brace->force_accept = 1;

        if (nbody_ode->length != ri_brace->encounter_N*3*2){
            // Just re-create the ODE
            printf("recreate ode\n");
            reb_ode_free(nbody_ode);
            nbody_ode = reb_ode_create(r, ri_brace->encounter_N*3*2);
            nbody_ode->derivatives = reb_integrator_brace_nbody_derivatives_barycentric;
            nbody_ode->needs_nbody = 0;
            r->ri_bs.first_or_last_step = 1;
        }
    }

    // Restore odes
    reb_ode_free(nbody_ode);
    free(r->odes);
    r->odes = odes_backup;
    r->N_allocated_odes = N_allocated_odes_backup;
    r->N_odes = N_odes_backup;

    r->t = old_t;

    // Resetting BS here reduces binary file size.
    reb_integrator_bs_reset(r);
}

void reb_integrator_brace_kepler_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    memcpy(ri_brace->particles_backup_kepler,r->particles,r->N*sizeof(struct reb_particle));
    reb_integrator_brace_whfast_step(r, _dt);
    reb_integrator_brace_bs_step(r, _dt);
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    particles[0].x  = 0; particles[0].y  = 0; particles[0].z  = 0;
    particles[0].vx = 0; particles[0].vy = 0; particles[0].vz = 0;
    for (int i=1;i<N;i++){
        const double f = particles[i].m/particles[0].m;
        particles[0].x  -= f * particles[i].x;
        particles[0].y  -= f * particles[i].y;
        particles[0].z  -= f * particles[i].z;
        particles[0].vx -= f * particles[i].vx;
        particles[0].vy -= f * particles[i].vy;
        particles[0].vz -= f * particles[i].vz;
    }
}


void reb_integrator_brace_part1(struct reb_simulation* r){
    // Do memory management and consistency checks in part1.
    // Actual integration is happening in part2.
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    const int N = r->N;
    
    if (r->N_var_config){
        reb_simulation_warning(r,"BRACE does not work with variational equations.");
    }

    if (ri_brace->N_allocated<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_brace->particles_backup       = realloc(ri_brace->particles_backup,sizeof(struct reb_particle)*N);
        ri_brace->particles_backup_kepler   = realloc(ri_brace->particles_backup_kepler,sizeof(struct reb_particle)*N);
        ri_brace->current_Ks             = realloc(ri_brace->current_Ks,sizeof(int)*N*N);
        ri_brace->current_Cs             = realloc(ri_brace->current_Cs,sizeof(int)*N);
        ri_brace->encounter_map          = realloc(ri_brace->encounter_map,sizeof(int)*N);
        ri_brace->N_allocated = N;
    }

    // Calculate collisions only with DIRECT or LINE method
    if (r->collision != REB_COLLISION_NONE && (r->collision != REB_COLLISION_DIRECT && r->collision != REB_COLLISION_LINE)){
        reb_simulation_warning(r,"BRACE only works with a direct or line collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_BRACE){
        reb_simulation_warning(r,"BRACE has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_BRACE;
    ri_brace->mode = REB_BRACE_MODE_NONE; // Do not calculate gravity in-between timesteps. BRACE will call reb_update_acceleration itself.

}

void reb_integrator_brace_pre_ts_check(struct reb_simulation* const r){
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;
    int (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = ri_brace->S ? ri_brace->S : reb_integrator_brace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = ri_brace->S_peri ? ri_brace->S_peri : reb_integrator_brace_switch_peri_default;
   
    // Clear encounter map
    for (unsigned int i=0; i<r->N; i++){
        ri_brace->encounter_map[i] = 0;
    }

    ri_brace->encounter_N = 0;

    for (int i=0; i<N; i++){
        ri_brace->current_Cs[i] = 0;
        for (unsigned int j = i + 1; j < N; j++){
            ri_brace->current_Ks[i*N+j] = 0;
        }
    }

    // Check for pericenter CE
    for (int i=1; i<Nactive; i++){
        if (_switch_peri(r, i)){
            ri_brace->current_Cs[i] = 1;
            printf("pre peri encounter\n");
        }
    }
    
    // Body-body
    for (int i=1; i<Nactive; i++){
        for (int j=i+1; j<N; j++){
            if (_switch(r, i, j)){
                ri_brace->current_Ks[i*N+j] = 1;
                if (ri_brace->encounter_map[i] == 0){
                    ri_brace->encounter_map[i] = 1; // trigger encounter
                    ri_brace->encounter_N++;
                }
                if (ri_brace->encounter_map[j] == 0){
                    ri_brace->encounter_map[j] = 1; // trigger encounter
                    ri_brace->encounter_N++;
                }
            }
        }
    }
}

int reb_integrator_brace_post_ts_check(struct reb_simulation* const r){
    // This function returns 1 if any new encounters occured.
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;
    int (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = ri_brace->S ? ri_brace->S : reb_integrator_brace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = ri_brace->S_peri ? ri_brace->S_peri : reb_integrator_brace_switch_peri_default;
    int new_close_encounter = 0;
    
    // Clear encounter maps
    for (unsigned int i=0; i<r->N; i++){
        ri_brace->encounter_map[i] = 0;
    }
    ri_brace->encounter_N = 0;
    
    // Check for pericenter CE if not already triggered from pre-timstep.
    for (int i=1; i<Nactive; i++){
        if (_switch_peri(r, i)){
            if (!ri_brace->current_Cs[i]){
                ri_brace->current_Cs[i] = 1;
                printf("new post peri encounter\n");
                new_close_encounter = 1;
            }
        }
    }

    // Body-body
    for (int i=1; i<Nactive; i++){
        for (int j=i+1; j<N; j++){
            if (_switch(r, i, j)){
                if (ri_brace->current_Ks[i*N+j] == 0){
                      new_close_encounter = 1;
                }
                ri_brace->current_Ks[i*N+j] = 1;
                if (ri_brace->encounter_map[i] == 0){
                    ri_brace->encounter_map[i] = 1; // trigger encounter
                    ri_brace->encounter_N++;
                }
                if (ri_brace->encounter_map[j] == 0){
                    ri_brace->encounter_map[j] = 1; // trigger encounter
                    ri_brace->encounter_N++;
                }
            }
        }
    }
    
    return new_close_encounter;
}

static void reb_integrator_brace_step(struct reb_simulation* const r){
    reb_integrator_brace_interaction_step(r, r->dt/2.);
    reb_integrator_brace_kepler_step(r, r->dt);
    reb_integrator_brace_com_step(r,r->dt);
    reb_integrator_brace_interaction_step(r, r->dt/2.);
}

void reb_integrator_brace_part2(struct reb_simulation* const r){
    struct reb_integrator_brace* const ri_brace = &(r->ri_brace);
    const int N = r->N;
    
    reb_integrator_brace_inertial_to_barycentric(r);
    
    // Create copy of all particle to allow for the step to be rejected.
    memcpy(ri_brace->particles_backup, r->particles, N*sizeof(struct reb_particle));
                        
    // This will be set to 1 if a collision occured.
    ri_brace->force_accept = 0;

    // Check if there are any close encounters
    reb_integrator_brace_pre_ts_check(r);
    
    // Attempt one step. 
    reb_integrator_brace_step(r);

    // We alaways accept the step if a collision occured as it is impossible to undo the collision.
    if (!ri_brace->force_accept){
        // We check again for close encounters to ensure time reversibility. 
        if (reb_integrator_brace_post_ts_check(r)){
            // New encounters were found. Will reject the step.
            // Revert particles to the beginning of the step.
            memcpy(r->particles, ri_brace->particles_backup, N*sizeof(struct reb_particle));

            // Do step again
            reb_integrator_brace_step(r);
        }
    }
    
    reb_integrator_brace_barycentric_to_inertial(r);
    
    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

void reb_integrator_brace_synchronize(struct reb_simulation* r){
}

void reb_integrator_brace_reset(struct reb_simulation* r){
    r->ri_brace.mode = REB_BRACE_MODE_NONE;
    r->ri_brace.encounter_N = 0;
    r->ri_brace.encounter_N_active = 0;
    r->ri_brace.r_crit_hill = 3;
    r->ri_brace.peri_crit_eta = 1.0;
    r->ri_brace.force_accept = 0;

    // Internal arrays (only used within one timestep)
    free(r->ri_brace.particles_backup);
    r->ri_brace.particles_backup = NULL;
    free(r->ri_brace.particles_backup_kepler);
    r->ri_brace.particles_backup_kepler = NULL;
    free(r->ri_brace.particles_backup_additional_forces);
    r->ri_brace.particles_backup_additional_forces = NULL;

    free(r->ri_brace.encounter_map);
    r->ri_brace.encounter_map = NULL;

    free(r->ri_brace.current_Ks);
    r->ri_brace.current_Ks = NULL;
    free(r->ri_brace.current_Cs);
    r->ri_brace.current_Cs = NULL;

    r->ri_brace.S = NULL;
    r->ri_brace.S_peri = NULL;
    
    r->ri_brace.peri_mode = REB_BRACE_PERI_FULL_BS;
    
    r->ri_brace.N_allocated = 0;
    r->ri_brace.N_allocated_additional_forces = 0;
}
