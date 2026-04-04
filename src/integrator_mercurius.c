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

#include "rebound.h"
#include "rebound_internal.h"
#include <string.h>
#include <math.h>
#include "gravity.h"
#include "integrator_mercurius.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

void reb_integrator_mercurius_calculate_acceleration_mode_encounter(struct reb_simulation* r);

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
    (void)r; // The simulation object is not used here.
}

double reb_integrator_mercurius_L_C4(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function C4 proposed by Hernandez (2019)
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return (70.*y*y*y*y -315.*y*y*y +540.*y*y -420.*y +126.)*y*y*y*y*y;
    }
    (void)r; // The simulation object is not used here.
}

double reb_integrator_mercurius_L_C5(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function C5 proposed by Hernandez (2019)
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return (-252.*y*y*y*y*y +1386.*y*y*y*y -3080.*y*y*y +3465.*y*y -1980.*y +462.)*y*y*y*y*y*y;
    }
    (void)r; // The simulation object is not used here.
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
    (void)r; // The simulation object is not used here.
}


void reb_integrator_mercurius_inertial_to_dh(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_vec3d com_pos = {0};
    struct reb_vec3d com_vel = {0};
    double mtot = 0.;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?r->N:r->N_active;
    const size_t N = r->N;
    for (size_t i=0;i<N_active;i++){
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
    struct reb_particle p0 = particles[0];
    for (size_t i=0;i<N;i++){ 
        particles[i].x -= p0.x;
        particles[i].y -= p0.y;
        particles[i].z -= p0.z;
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
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?r->N:r->N_active;
    for (size_t i=1;i<N_active;i++){
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
    // Should not matter unless collisions occurred.
    particles[0].x = r->ri_mercurius.com_pos.x - temp.x; 
    particles[0].y = r->ri_mercurius.com_pos.y - temp.y; 
    particles[0].z = r->ri_mercurius.com_pos.z - temp.z; 

    for (size_t i=1;i<N;i++){
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
    struct reb_integrator_mercurius* rim = &(r->ri_mercurius);
    struct reb_particle* const particles = r->particles;
    struct reb_particle* const particles_backup = rim->particles_backup;
    const double* const dcrit = rim->dcrit;
    const size_t N = r->N;
    const size_t N_active = r->N_active==SIZE_MAX?r->N:r->N_active;
    const double dt = r->dt;
    rim->encounter_N = 1;
    rim->encounter_map[0] = 1;
    if (r->testparticle_type==1){
        rim->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        rim->tponly_encounter = 1;
    }
    for (size_t i=1; i<N; i++){
        rim->encounter_map[i] = 0;
    }
    for (size_t i=0; i<N_active; i++){
        for (size_t j=i+1; j<N; j++){
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

            double dcritmax2 = MAX(dcrit[i],dcrit[j]);
            dcritmax2 *= 1.21*dcritmax2;
            if (rmin < dcritmax2){
                if (rim->encounter_map[i]==0){
                    rim->encounter_map[i] = i;
                    rim->encounter_N++;
                }
                if (rim->encounter_map[j]==0){
                    rim->encounter_map[j] = j;
                    rim->encounter_N++;
                }
                if (j<N_active){ // Two massive particles have a close encounter
                    rim->tponly_encounter = 0;
                }
            }
        }
    }
}

void reb_integrator_mercurius_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    for (size_t i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_mercurius_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N_active = r->N_active==SIZE_MAX?r->N:r->N_active;
    const size_t N = r->testparticle_type==0 ? N_active: r->N;
    double px=0., py=0., pz=0.;
    for (size_t i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m; 
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px /= r->particles[0].m;
    py /= r->particles[0].m;
    pz /= r->particles[0].m;
    const size_t N_all = r->N;
    for (size_t i=1;i<N_all;i++){
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
    const size_t N = r->N;
    int timestep_too_large = 0;
    for (size_t i=1;i<N;i++){
        timestep_too_large |= reb_integrator_whfast_kepler_solver(&particles[i],r->G*particles[0].m,dt,NULL); // in dh
    }
    if (timestep_too_large && r->ri_whfast.timestep_warning == 0){
        r->ri_whfast.timestep_warning++;
        reb_simulation_warning(r,"Kepler solver convergence issue. Timestep is larger than at least one orbital period.");
    }
}

static void reb_mercurius_encounter_step(struct reb_simulation* const r, const double _dt){
    // Only particles having a close encounter are integrated by IAS15.
    struct reb_integrator_mercurius* rim = &(r->ri_mercurius);
    if (rim->encounter_N<2){
        return; // If there are no particles (other than the star) having a close encounter, then there is nothing to do.
    }
    size_t N_active = (r->N_active==SIZE_MAX)?r->N:r->N_active;
    size_t i_enc = 0;
    rim->encounter_N_active = 0;
    for (size_t i=0; i<r->N; i++){
        if(rim->encounter_map[i]){  
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = rim->particles_backup[i];     // Use coordinates before whfast step
            rim->encounter_map[i_enc] = i;
            i_enc++;
            if (i<N_active){
                rim->encounter_N_active++;
                if (rim->tponly_encounter){
                    rim->particles_backup[i] = tmp;         // Make copy of particles after the kepler step.
                                                            // used to restore the massive objects' states in the case
                                                            // of only massless test-particle encounters
                }
            }
        }
    }

    rim->mode = REB_MERCURIUS_MODE_ENCOUNTER;

    // run
    const double old_dt = r->dt;
    const double dtsign = old_dt>=0.?1.:-1.;
    const double old_t = r->t;
    double t_needed = r->t + _dt; 

    reb_integrator_ias15_reset(r);

    r->dt = 0.0001*_dt; // start with a small timestep.

    // No additional forces during encounter
    r->gravity = REB_GRAVITY_CUSTOM;
    r->gravity_custom = reb_integrator_mercurius_calculate_acceleration_mode_encounter;

    while(dtsign*r->t < dtsign*t_needed && fabs(r->dt/old_dt)>1e-14 && r->status<=0){
        struct reb_particle star = r->particles[0]; // backup velocity
        r->particles[0].vx = 0; // star does not move in dh 
        r->particles[0].vy = 0;
        r->particles[0].vz = 0;
        reb_integrator_ias15_step(r);
        r->particles[0].vx = star.vx; // restore every timestep for collisions
        r->particles[0].vy = star.vy;
        r->particles[0].vz = star.vz;

        if (dtsign*(r->t+r->dt) > dtsign*t_needed){
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

        struct reb_particle p0 = r->particles[0];
        star.vx = p0.vx; // keep track of changed star velocity for later collisions
        star.vy = p0.vy;
        star.vz = p0.vz;
        if (p0.x !=0 || p0.y !=0 || p0.z !=0){
            // Collision with star occurred
            // Shift all particles back to heliocentric coordinates
            // Ignore stars velocity:
            //   - will not be used after this
            //   - com velocity is unchanged. this velocity will be used
            //     to reconstruct star's velocity later.
            for (size_t i=0; i<r->N; i++){
                r->particles[i].x -= p0.x;
                r->particles[i].y -= p0.y;
                r->particles[i].z -= p0.z;
            }
        }
    }

    // if only test particles encountered massive bodies, reset the
    // massive body coordinates to their post Kepler step state
    if(rim->tponly_encounter){
        for (size_t i=1;i<rim->encounter_N_active;i++){
            size_t mi = rim->encounter_map[i];
            r->particles[mi] = rim->particles_backup[mi];
        }
    }

    // Reset constant for global particles
    r->t = old_t;
    r->dt = old_dt;
    rim->mode = REB_MERCURIUS_MODE_WH;

}

double reb_integrator_mercurius_calculate_dcrit_for_particle(struct reb_simulation* r, size_t i){
    struct reb_integrator_mercurius* const rim = &(r->ri_mercurius);
    const double m0 = r->particles[0].m;
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
    dcrit = MAX(dcrit, rim->r_crit_hill*a*cbrt(r->particles[i].m/(3.*r->particles[0].m)));
    // Criteria 4: physical radius
    dcrit = MAX(dcrit, 2.*r->particles[i].r);
    return dcrit;
}
        
void reb_integrator_mercurius_calculate_acceleration_mode_encounter(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = r->ri_mercurius.L;
    const double m0 = r->particles[0].m;
    const double* const dcrit = r->ri_mercurius.dcrit;
    const size_t encounter_N = r->ri_mercurius.encounter_N;
    const size_t encounter_N_active = r->ri_mercurius.encounter_N_active;
    size_t* map = r->ri_mercurius.encounter_map;
#ifndef OPENMP
    particles[0].ax = 0; // map[0] is always 0 
    particles[0].ay = 0; 
    particles[0].az = 0; 
    // Acceleration due to star
    for (size_t i=1; i<encounter_N; i++){
        size_t mi = map[i];
        const double x = particles[mi].x;
        const double y = particles[mi].y;
        const double z = particles[mi].z;
        const double _r = sqrt(x*x + y*y + z*z + softening2);
        double prefact = -G/(_r*_r*_r)*m0;
        particles[mi].ax    = prefact*x;
        particles[mi].ay    = prefact*y;
        particles[mi].az    = prefact*z;
    }
    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
    // Interactions between active-active
    for (size_t i=2; i<encounter_N_active; i++){
        size_t mi = map[i];
        for (size_t j=1; j<i; j++){
            size_t mj = map[j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double L = _L(r,_r,dcritmax);
            double prefact = G*(1.-L)/(_r*_r*_r);
            double prefactj = -prefact*particles[mj].m;
            double prefacti = prefact*particles[mi].m;
            particles[mi].ax    += prefactj*dx;
            particles[mi].ay    += prefactj*dy;
            particles[mi].az    += prefactj*dz;
            particles[mj].ax    += prefacti*dx;
            particles[mj].ay    += prefacti*dy;
            particles[mj].az    += prefacti*dz;
        }
    }
    // Interactions between active-testparticle
    const size_t startitestp = MAX(encounter_N_active,2);
    for (size_t i=startitestp; i<encounter_N; i++){
        size_t mi = map[i];
        for (size_t j=1; j<encounter_N_active; j++){
            size_t mj = map[j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double L = _L(r,_r,dcritmax);
            double prefact = G*(1.-L)/(_r*_r*_r);
            double prefactj = -prefact*particles[mj].m;
            particles[mi].ax    += prefactj*dx;
            particles[mi].ay    += prefactj*dy;
            particles[mi].az    += prefactj*dz;
            if (_testparticle_type){
                double prefacti = prefact*particles[mi].m;
                particles[mj].ax    += prefacti*dx;
                particles[mj].ay    += prefacti*dy;
                particles[mj].az    += prefacti*dz;
            }
        }
    }
#else // OPENMP
    particles[0].ax = 0; // map[0] is always 0 
    particles[0].ay = 0; 
    particles[0].az = 0; 
    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
#pragma omp parallel for schedule(guided)
    for (size_t i=1; i<encounter_N; i++){
        size_t mi = map[i];
        particles[mi].ax = 0; 
        particles[mi].ay = 0; 
        particles[mi].az = 0; 
        // Acceleration due to star
        const double x = particles[mi].x;
        const double y = particles[mi].y;
        const double z = particles[mi].z;
        const double _r = sqrt(x*x + y*y + z*z + softening2);
        double prefact = -G/(_r*_r*_r)*m0;
        particles[mi].ax    += prefact*x;
        particles[mi].ay    += prefact*y;
        particles[mi].az    += prefact*z;
        for (size_t j=1; j<encounter_N_active; j++){
            if (i==j) continue;
            size_t mj = map[j];
            const double dx = x - particles[mj].x;
            const double dy = y - particles[mj].y;
            const double dz = z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double L = _L(r,_r,dcritmax);
            double prefact = -G*particles[mj].m*(1.-L)/(_r*_r*_r);
            particles[mi].ax    += prefact*dx;
            particles[mi].ay    += prefact*dy;
            particles[mi].az    += prefact*dz;
        }
    }
    if (_testparticle_type){
#pragma omp parallel for schedule(guided)
        for (size_t i=1; i<encounter_N_active; i++){
            size_t mi = map[i];
            const double x = particles[mi].x;
            const double y = particles[mi].y;
            const double z = particles[mi].z;
            for (size_t j=encounter_N_active; j<encounter_N; j++){
                size_t mj = map[j];
                const double dx = x - particles[mj].x;
                const double dy = y - particles[mj].y;
                const double dz = z - particles[mj].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                const double L = _L(r,_r,dcritmax);
                double prefact = -G*particles[mj].m*(1.-L)/(_r*_r*_r);
                particles[mi].ax    += prefact*dx;
                particles[mi].ay    += prefact*dy;
                particles[mi].az    += prefact*dz;
            }
        }
    }
#endif // OPENMP
}

static void reb_integrator_mercurius_calculate_acceleration_mode_wh(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = r->ri_mercurius.L;
    const double* const dcrit = r->ri_mercurius.dcrit;
#ifndef OPENMP
    for (size_t i=0; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    for (size_t i=2; i<N_active; i++){
        if (reb_sigint > 1) return;
        for (size_t j=1; j<i; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            const double prefacti = prefact*particles[i].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            particles[j].ax    += prefacti*dx;
            particles[j].ay    += prefacti*dy;
            particles[j].az    += prefacti*dz;
        }
    }
    const size_t startitestp = MAX(N_active,2);
    for (size_t i=startitestp; i<N; i++){
        if (reb_sigint > 1) return;
        for (size_t j=1; j<N_active; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            if (_testparticle_type){
                const double prefacti = prefact*particles[i].m;
                particles[j].ax    += prefacti*dx;
                particles[j].ay    += prefacti*dy;
                particles[j].az    += prefacti*dz;
            }
        }
    }
#else // OPENMP
    particles[0].ax = 0; 
    particles[0].ay = 0; 
    particles[0].az = 0; 
#pragma omp parallel for schedule(guided)
    for (size_t i=1; i<N; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
        for (size_t j=1; j<N_active; j++){
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,_r,dcritmax);
            const double prefact = -G*particles[j].m*L/(_r*_r*_r);
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    if (_testparticle_type){
        for (size_t i=1; i<N_active; i++){
            for (size_t j=N_active; j<N; j++){
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                const double dcritmax = MAX(dcrit[i],dcrit[j]);
                const double L = _L(r,_r,dcritmax);
                const double prefact = -G*particles[j].m*L/(_r*_r*_r);
                particles[i].ax    += prefact*dx;
                particles[i].ay    += prefact*dy;
                particles[i].az    += prefact*dz;
            }
        }
    }
#endif // OPENMP

    if (r->additional_forces){
        // Additional forces are only calculated in the kick step, not during close encounter
        // shift pos and velocity so that external forces are calculated in inertial frame
        // Note: Copying avoids degrading floating point performance
        if(r->N>r->ri_mercurius.N_allocated_additional_forces){
            r->ri_mercurius.particles_backup_additional_forces = realloc(r->ri_mercurius.particles_backup_additional_forces, r->N*sizeof(struct reb_particle));
            r->ri_mercurius.N_allocated_additional_forces = r->N;
        }
        memcpy(r->ri_mercurius.particles_backup_additional_forces,r->particles,r->N*sizeof(struct reb_particle)); 
        reb_integrator_mercurius_dh_to_inertial(r);

        r->additional_forces(r);

        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle* restrict const backup = r->ri_mercurius.particles_backup_additional_forces;
        for (size_t i=0;i<r->N;i++){
            particles[i].x = backup[i].x;
            particles[i].y = backup[i].y;
            particles[i].z = backup[i].z;
            particles[i].vx = backup[i].vx;
            particles[i].vy = backup[i].vy;
            particles[i].vz = backup[i].vz;
        }
    }
}


void reb_integrator_mercurius_step(struct reb_simulation* r){
    if (r->N_var_config){
        reb_simulation_warning(r,"Mercurius does not work with variational equations.");
    }

    struct reb_integrator_mercurius* const rim = &(r->ri_mercurius);
    const size_t N = r->N;

    if (rim->N_allocated_dcrit<N){
        // Need to safe these arrays in Simulationarchive
        rim->dcrit              = realloc(rim->dcrit, sizeof(double)*N);
        rim->N_allocated_dcrit = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_r_crit_this_timestep        = 1;
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
        rim->recalculate_coordinates_this_timestep = 1;
    }
    if (rim->N_allocated<N){
        // These arrays are only used within one timestep. 
        // Can be recreated without loosing bit-wise reproducibility
        rim->particles_backup   = realloc(rim->particles_backup,sizeof(struct reb_particle)*N);
        rim->encounter_map      = realloc(rim->encounter_map,sizeof(size_t)*N);
        rim->N_allocated = N;
    }
    if (rim->safe_mode || rim->recalculate_coordinates_this_timestep){
        if (rim->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r);
            reb_simulation_warning(r,"MERCURIUS: Recalculating heliocentric coordinates but coordinates were not synchronized before.");
        }
        reb_integrator_mercurius_inertial_to_dh(r);
        rim->recalculate_coordinates_this_timestep = 0;
    }

    if (rim->recalculate_r_crit_this_timestep){
        rim->recalculate_r_crit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r);
            reb_integrator_mercurius_inertial_to_dh(r);
            rim->recalculate_coordinates_this_timestep = 0;
            reb_simulation_warning(r,"MERCURIUS: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        rim->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        for (size_t i=1;i<N;i++){
            rim->dcrit[i] = reb_integrator_mercurius_calculate_dcrit_for_particle(r, i);
        }
    }

    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_simulation_warning(r,"Mercurius only works with a direct collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_CUSTOM){
        reb_simulation_warning(r,"Mercurius has its own gravity routine. Gravity routine set by the user will be ignored.");
    }
    rim->mode = REB_MERCURIUS_MODE_WH;

    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurius_L_mercury;
    }

    reb_integrator_mercurius_calculate_acceleration_mode_wh(r);

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
    struct reb_integrator_mercurius* const rim = &(r->ri_mercurius);
    if (rim->is_synchronized == 0){
        rim->mode = REB_MERCURIUS_MODE_WH;
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurius_L_mercury;
        }
        reb_integrator_mercurius_calculate_acceleration_mode_wh(r);
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);

        reb_integrator_mercurius_dh_to_inertial(r);

        rim->recalculate_coordinates_this_timestep = 1; 
        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurius_reset(struct reb_simulation* r){
    r->ri_mercurius.L = NULL;
    r->ri_mercurius.mode = REB_MERCURIUS_MODE_WH;
    r->ri_mercurius.encounter_N = 0;
    r->ri_mercurius.encounter_N_active = 0;
    r->ri_mercurius.r_crit_hill = 3;
    r->ri_mercurius.tponly_encounter = 0;
    r->ri_mercurius.safe_mode = 1;
    r->ri_mercurius.is_synchronized = 1;
    r->ri_mercurius.recalculate_coordinates_this_timestep = 0;
    r->ri_mercurius.recalculate_r_crit_this_timestep = 0;
    // Internal arrays (only used within one timestep)
    free(r->ri_mercurius.particles_backup);
    r->ri_mercurius.particles_backup = NULL;
    free(r->ri_mercurius.particles_backup_additional_forces);
    r->ri_mercurius.particles_backup_additional_forces = NULL;
    free(r->ri_mercurius.encounter_map);
    r->ri_mercurius.encounter_map = NULL;
    r->ri_mercurius.N_allocated = 0;
    r->ri_mercurius.N_allocated_additional_forces = 0;
    // dcrit array
    free(r->ri_mercurius.dcrit);
    r->ri_mercurius.dcrit = NULL;
    r->ri_mercurius.N_allocated_dcrit = 0;
}

