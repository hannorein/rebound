/**
 * @file    integrator_trace.c
 * @brief   TRACE
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

#include "rebound.h"
#include "rebound_internal.h"
#include <time.h>
#include <string.h>
#include <math.h>
#include "gravity.h"
#include "integrator_trace.h"
#include "integrator_whfast.h"
#include "integrator_bs.h"
#include "integrator_ias15.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

const struct reb_integrator reb_integrator_trace = {
    .id = 25,
    .step = reb_integrator_trace_step,
    .reset = reb_integrator_trace_reset,
    .will_remove_particle = reb_integrator_trace_will_remove_particle,
    .did_add_particle = reb_integrator_trace_did_add_particle,
};

int reb_integrator_trace_switch_default(struct reb_simulation* const r, const size_t i, const size_t j){
    // Returns 1 for close encounter between i and j, 0 otherwise
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
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

    double r_crit_hill2 = ri_trace->r_crit_hill*ri_trace->r_crit_hill;
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

int reb_integrator_trace_switch_peri_default(struct reb_simulation* const r, const size_t j){
    // Following Pham et al (2024)
    const struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
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
    double dt_prs2 = ri_trace->peri_crit_eta * ri_trace->peri_crit_eta * tau_prs2;

    if (r->dt * r->dt > dt_prs2){
        return 1;
    }else{
        return 0;
    }
}

int reb_integrator_trace_switch_peri_none(struct reb_simulation* const r, const size_t j){
    // No pericenter flags
    (void)r; // Not used
    (void)j; // Not used
    return 0;
}

void reb_integrator_trace_inertial_to_dh(struct reb_simulation* r){
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
    r->ri_trace.com_pos = com_pos;
    r->ri_trace.com_vel = com_vel;
}

void reb_integrator_trace_dh_to_inertial(struct reb_simulation* r){
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
    particles[0].x = r->ri_trace.com_pos.x - temp.x;
    particles[0].y = r->ri_trace.com_pos.y - temp.y;
    particles[0].z = r->ri_trace.com_pos.z - temp.z;

    for (size_t i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += r->ri_trace.com_vel.x;
        particles[i].vy += r->ri_trace.com_vel.y;
        particles[i].vz += r->ri_trace.com_vel.z;
    }
    particles[0].vx = r->ri_trace.com_vel.x - temp.vx;
    particles[0].vy = r->ri_trace.com_vel.y - temp.vy;
    particles[0].vz = r->ri_trace.com_vel.z - temp.vz;
}

static void reb_integrator_trace_calculate_acceleration_mode_interaction(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const int _testparticle_type   = r->testparticle_type;
#ifndef OPENMP
    for (size_t i=0; i<N; i++){
        particles[i].ax = 0;
        particles[i].ay = 0;
        particles[i].az = 0;
    }
    for (size_t i=2; i<N_active; i++){
        if (reb_sigint > 1) return;
        for (size_t j=1; j<i; j++){
            if (r->ri_trace.current_Ks[j*N+i]) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double prefact = G / (_r*_r*_r);
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
            if (r->ri_trace.current_Ks[j*N+i]) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double prefact = G / (_r*_r*_r);
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
            if (r->ri_trace.current_Ks[j*N+i]) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double prefact = -G*particles[j].m/(_r*_r*_r);
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    if (_testparticle_type){
        for (size_t i=1; i<N_active; i++){
            for (size_t j=N_active; j<N; j++){
                if (r->ri_trace.current_Ks[j*N+i]) continue;
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                const double prefact = -G*particles[j].m/(_r*_r*_r);
                particles[i].ax    += prefact*dx;
                particles[i].ay    += prefact*dy;
                particles[i].az    += prefact*dz;
            }
        }
    }
#endif // OPENMP

    // Handle Additional forces
    if (r->additional_forces){
        // shift pos and velocity so that external forces are calculated in inertial frame
        // Note: Copying avoids degrading floating point performance
        // We should NOT do this in FULL mode, already in inertial frame
        if(r->N>r->ri_trace.N_allocated_additional_forces){
            r->ri_trace.particles_backup_additional_forces = realloc(r->ri_trace.particles_backup_additional_forces, r->N*sizeof(struct reb_particle));
            r->ri_trace.N_allocated_additional_forces = r->N;
        }
        memcpy(r->ri_trace.particles_backup_additional_forces,r->particles,r->N*sizeof(struct reb_particle));
        reb_integrator_trace_dh_to_inertial(r);
        r->additional_forces(r);
        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle* restrict const backup = r->ri_trace.particles_backup_additional_forces;
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

static void reb_integrator_trace_calculate_acceleration_mode_kepler(struct reb_simulation* r){
    // Kepler Step
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    const double m0 = r->particles[0].m;
    const size_t encounter_N = r->ri_trace.encounter_N;
    const size_t encounter_N_active = r->ri_trace.encounter_N_active;
    size_t* map = r->ri_trace.encounter_map;
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
        double prefact = -G * m0 / (_r*_r*_r);
        particles[mi].ax    = prefact*x;
        particles[mi].ay    = prefact*y;
        particles[mi].az    = prefact*z;
    }

    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
    // Interactions between active-active
    if (encounter_N_active > 2){ // if two or less, no active-active planets
        for (size_t i=2; i<encounter_N_active; i++){
            size_t mi = map[i];
            for (size_t j=1; j<i; j++){
                size_t mj = map[j];
                if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
                const double dx = particles[mi].x - particles[mj].x;
                const double dy = particles[mi].y - particles[mj].y;
                const double dz = particles[mi].z - particles[mj].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                double prefact = G/(_r*_r*_r);
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
    }

    // Interactions between active-testparticle
    const size_t startitestp = MAX(encounter_N_active,2);
    for (size_t i=startitestp; i<encounter_N; i++){
        size_t mi = map[i];
        for (size_t j=1; j<encounter_N_active; j++){
            size_t mj = map[j];
            if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            double prefact = G/(_r*_r*_r);
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
            if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
            const double dx = x - particles[mj].x;
            const double dy = y - particles[mj].y;
            const double dz = z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            double prefact = -G*particles[mj].m/(_r*_r*_r);
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
                if (!r->ri_trace.current_Ks[mj*N+mi]) continue;
                const double dx = x - particles[mj].x;
                const double dy = y - particles[mj].y;
                const double dz = z - particles[mj].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                double prefact = -G*particles[mj].m/(_r*_r*_r);
                particles[mi].ax    += prefact*dx;
                particles[mi].ay    += prefact*dy;
                particles[mi].az    += prefact*dz;
            }
        }
    }
#endif // OPENMP
}

void reb_integrator_trace_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    r->ri_trace.mode = REB_TRACE_MODE_INTERACTION;
    reb_integrator_trace_calculate_acceleration_mode_interaction(r);
    for (size_t i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_trace_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;

    struct reb_integrator_trace* ri_trace = &(r->ri_trace);
    const int current_C = ri_trace->current_C;
    if (current_C) return; // No jump step for pericenter approaches

    const size_t N_active = r->N_active==SIZE_MAX?r->N:r->N_active;

    // If TP type 1, use r->N. Else, use N_active.
    const size_t N = r->testparticle_type==0 ? N_active: r->N;

    double px=0., py=0., pz=0.;
    for (size_t i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m;
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px *= dt/r->particles[0].m;
    py *= dt/r->particles[0].m;
    pz *= dt/r->particles[0].m;

    const size_t N_all = r->N;
    for (size_t i=1;i<N_all;i++){
        particles[i].x += px;
        particles[i].y += py;
        particles[i].z += pz;
    }
}

void reb_integrator_trace_com_step(struct reb_simulation* const r, double dt){
    r->ri_trace.com_pos.x += dt*r->ri_trace.com_vel.x;
    r->ri_trace.com_pos.y += dt*r->ri_trace.com_vel.y;
    r->ri_trace.com_pos.z += dt*r->ri_trace.com_vel.z;
}

void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const size_t N = r->N;
    int timestep_too_large = 0;
    for (size_t i=1;i<N;i++){
        timestep_too_large |= reb_integrator_whfast_kepler_solver(&particles[i],r->G*particles[0].m,dt,NULL);
    }
    if (timestep_too_large && r->ri_whfast.timestep_warning == 0){
        r->ri_whfast.timestep_warning++;
        reb_simulation_warning(r,"Kepler solver convergence issue. Timestep is larger than at least one orbital period.");
    }
}

void reb_integrator_trace_update_particles(struct reb_simulation* r, const double* y){
    size_t N = r->ri_trace.encounter_N;
    size_t* map = r->ri_trace.encounter_map;

    for (size_t i=0; i<N; i++){
        size_t mi = map[i];
        struct reb_particle* const p = &(r->particles[mi]);
        p->x  = y[i*6+0];
        p->y  = y[i*6+1];
        p->z  = y[i*6+2];
        p->vx = y[i*6+3];
        p->vy = y[i*6+4];
        p->vz = y[i*6+5];
    }
}

void reb_integrator_trace_nbody_derivatives(struct reb_ode* ode, double* const yDot, const double* const y, double const t){
    (void)t; // Not timedependent.
    struct reb_simulation* const r = ode->r;
    // TRACE always needs this to ensure the right Hamiltonian is evolved
    reb_integrator_trace_update_particles(r, y);
    reb_integrator_trace_calculate_acceleration_mode_kepler(r);

    double px=0., py=0., pz=0.;
    size_t* map = r->ri_trace.encounter_map;
    size_t N = r->ri_trace.encounter_N;

    if (map==NULL){
        reb_simulation_error(r, "Cannot access TRACE map from BS.");
        return;
    }

    // Kepler Step
    // This is only for pericenter approach
    if (r->ri_trace.current_C){
        for (size_t i=1;i<r->N;i++){ // all particles
            px += r->particles[i].vx*r->particles[i].m; // in dh
            py += r->particles[i].vy*r->particles[i].m;
            pz += r->particles[i].vz*r->particles[i].m;
        }
        px /= r->particles[0].m;
        py /= r->particles[0].m;
        pz /= r->particles[0].m;

    }
    yDot[0*6+0] = 0.0;
    yDot[0*6+1] = 0.0;
    yDot[0*6+2] = 0.0;
    yDot[0*6+3] = 0.0;
    yDot[0*6+4] = 0.0;
    yDot[0*6+5] = 0.0;

    for (size_t i=1; i<N; i++){
        size_t mi = map[i];
        const struct reb_particle p = r->particles[mi];
        yDot[i*6+0] = p.vx + px; // Already checked for current_L
        yDot[i*6+1] = p.vy + py;
        yDot[i*6+2] = p.vz + pz;
        yDot[i*6+3] = p.ax;
        yDot[i*6+4] = p.ay;
        yDot[i*6+5] = p.az;
    }
}

void reb_integrator_trace_bs_step(struct reb_simulation* const r, double dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);

    if (ri_trace->encounter_N < 2){
        // No close encounters, skip
        return;
    }

    size_t i_enc = 0;
    const size_t N_active = r->N_active==SIZE_MAX ? r->N : r->N_active;
    ri_trace->encounter_N_active = 0;
    for (size_t i=0; i<r->N; i++){
        if(ri_trace->encounter_map[i]){
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = ri_trace->particles_backup_kepler[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            ri_trace->encounter_map[i_enc] = i;
            i_enc++;
            if (i<N_active){
                ri_trace->encounter_N_active++;
                if (ri_trace->tponly_encounter){
                    ri_trace->particles_backup_kepler[i] = tmp;         // Make copy of particles after the kepler step.
                                                                        // used to restore the massive objects' states in the case
                                                                        // of only massless test-particle encounters
                }
            }
        }
    }

    ri_trace->mode = REB_TRACE_MODE_KEPLER;
    r->map = ri_trace->encounter_map; // for collision search
    r->N_map = ri_trace->encounter_N;
    r->gravity = REB_GRAVITY_CUSTOM;
    r->gravity_custom = reb_integrator_trace_calculate_acceleration_mode_kepler;

    // Only Partial BS uses this step 
    if (ri_trace->peri_mode == REB_TRACE_PERI_PARTIAL_BS || !ri_trace->current_C){
        // run
        const double old_dt = r->dt;
        const double old_t = r->t;
        const double t_needed = r->t + dt;
        reb_integrator_bs_reset(r);

        // Temporarily remove all odes for BS step
        struct reb_ode** odes_backup = r->odes;
        size_t N_allocated_odes_backup = r->N_allocated_odes;
        size_t N_odes_backup = r->N_odes;
        r->odes = NULL;
        r->N_allocated_odes = 0;
        r->N_odes = 0;

        // Temporarily add new nbody ode for BS step
        struct reb_ode* nbody_ode = reb_ode_create(r, ri_trace->encounter_N*3*2);
        nbody_ode->derivatives = reb_integrator_trace_nbody_derivatives;
        nbody_ode->needs_nbody = 0;

        // TODO: Support backwards integrations
        while(r->t < t_needed && fabs(dt/old_dt)>1e-14 && r->status<=0){
            double* y = nbody_ode->y;

            // In case of overshoot
            if (r->t + dt >  t_needed){
                dt = t_needed - r->t;
            }

            struct reb_particle star = r->particles[0]; // backup velocity
            r->particles[0].vx = 0; // star does not move in dh
            r->particles[0].vy = 0;
            r->particles[0].vz = 0;

            for (size_t i=0; i<ri_trace->encounter_N; i++){
                const size_t mi = ri_trace->encounter_map[i];
                const struct reb_particle p = r->particles[mi];
                y[i*6+0] = p.x;
                y[i*6+1] = p.y;
                y[i*6+2] = p.z;
                y[i*6+3] = p.vx;
                y[i*6+4] = p.vy;
                y[i*6+5] = p.vz;
            }

            int success = reb_integrator_bs_step_odes(r, dt);
            if (success){
                r->t += dt;
            }
            dt = r->ri_bs.dt_proposed;
            reb_integrator_trace_update_particles(r, nbody_ode->y);

            r->particles[0].vx = star.vx; // restore every timestep for collisions
            r->particles[0].vy = star.vy;
            r->particles[0].vz = star.vz;

            if (success){
                // Only do a collision search for accepted steps.
                reb_collision_search(r);
                if (r->N_collisions) r->ri_trace.force_accept = 1;
            }

            if (nbody_ode->length != ri_trace->encounter_N*3*2){
                // Just re-create the ODE
                reb_ode_free(nbody_ode);
                nbody_ode = reb_ode_create(r, ri_trace->encounter_N*3*2);
                nbody_ode->derivatives = reb_integrator_trace_nbody_derivatives;
                nbody_ode->needs_nbody = 0;
                r->ri_bs.first_or_last_step = 1;
            }

            struct reb_particle p0 = r->particles[0];
            star.vx = p0.vx; // keep track of changed star velocity for later collisions
            star.vy = p0.vy;
            star.vz = p0.vz;

            if (r->particles[0].x !=0 || r->particles[0].y !=0 || r->particles[0].z !=0){
                // Collision with star occurred
                // Shift all particles back to heliocentric coordinates
                // Ignore stars velocity:
                //   - will not be used after this
                //   - com velocity is unchanged. this velocity will be used
                //     to reconstruct star's velocity later.
                for (size_t i=0; i<r->N;i++){
                    r->particles[i].x -= p0.x;
                    r->particles[i].y -= p0.y;
                    r->particles[i].z -= p0.z;
                }
            }
        }

        // if only test particles encountered massive bodies, reset the
        // massive body coordinates to their post Kepler step state
        if(ri_trace->tponly_encounter){
            for (size_t i=1; i < ri_trace->encounter_N_active; i++){
                size_t mi = ri_trace->encounter_map[i];
                r->particles[mi] = ri_trace->particles_backup_kepler[mi];
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
    r->map = NULL; // for collision search
    r->N_map = 0;
}

void reb_integrator_trace_kepler_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    memcpy(ri_trace->particles_backup_kepler,r->particles,r->N*sizeof(struct reb_particle));
    reb_integrator_trace_whfast_step(r, _dt);
    reb_integrator_trace_bs_step(r, _dt);
}


void reb_integrator_trace_pre_ts_check(struct reb_simulation* const r){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const size_t N = r->N;
    const size_t Nactive = r->N_active==SIZE_MAX?r->N:r->N_active;
    int (*_switch) (struct reb_simulation* const r, const size_t i, const size_t j) = ri_trace->S ? ri_trace->S : reb_integrator_trace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const size_t j) = ri_trace->S_peri ? ri_trace->S_peri : reb_integrator_trace_switch_peri_default;

    // Clear encounter map
    for (size_t i=1; i<r->N; i++){
        ri_trace->encounter_map[i] = 0;
    }
    ri_trace->encounter_map[0] = 1;
    ri_trace->encounter_N = 1;

    // Reset encounter triggers.
    ri_trace->current_C = 0;

    for (size_t i = 0; i < N; i++){
        for (size_t j = i + 1; j < N; j++){
            ri_trace->current_Ks[i*N+j] = 0;
        }
    }

    if (r->testparticle_type == 1){
        ri_trace->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        ri_trace->tponly_encounter = 1;
    }

    // Check for pericenter CE
    for (size_t j = 1; j < Nactive; j++){
        if (_switch_peri(r, j)){
            ri_trace->current_C = 1;
            if (ri_trace->peri_mode == REB_TRACE_PERI_FULL_BS || ri_trace->peri_mode == REB_TRACE_PERI_FULL_IAS15){
                // Everything will be integrated with BS/IAS15. No need to check any further.
                return;
            }
            if (j < Nactive){ // Two massive particles have a close encounter
                ri_trace->tponly_encounter = 0;
                break; // No need to check other particles
            }
        }
    }

    if (ri_trace->current_C){
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        ri_trace->encounter_N = N;
        for (size_t i = 1; i < N; i++){
            ri_trace->encounter_map[i] = 1; //  trigger encounter
        }

    }

    // Body-body
    // there cannot be TP-TP CEs
    for (size_t i = 0; i < Nactive; i++){ // Check central body, for collisions
        for (size_t j = i + 1; j < N; j++){
            if (_switch(r, i, j)){
                ri_trace->current_Ks[i*N+j] = 1;
                if (ri_trace->encounter_map[i] == 0){
                    ri_trace->encounter_map[i] = 1; // trigger encounter
                    ri_trace->encounter_N++;
                }
                if (ri_trace->encounter_map[j] == 0){
                    ri_trace->encounter_map[j] = 1; // trigger encounter
                    ri_trace->encounter_N++;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    ri_trace->tponly_encounter = 0;
                }
            }
        }
    }
    memcpy(ri_trace->encounter_map_backup, ri_trace->encounter_map, N*sizeof(size_t));
}

double reb_integrator_trace_post_ts_check(struct reb_simulation* const r){
    // This function returns 1 if any new encounters occurred.
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const size_t N = r->N;
    const size_t Nactive = r->N_active==SIZE_MAX?r->N:r->N_active;
    int (*_switch) (struct reb_simulation* const r, const size_t i, const size_t j) = ri_trace->S ? ri_trace->S : reb_integrator_trace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const size_t j) = ri_trace->S_peri ? ri_trace->S_peri : reb_integrator_trace_switch_peri_default;
    size_t new_close_encounter = 0; // New CEs

    // Set this from pre-ts encounter map. I don't think we need to reset encounter_N here.
    memcpy(ri_trace->encounter_map, ri_trace->encounter_map_backup, N*sizeof(size_t));

    if (!ri_trace->current_C){
        // Check for pericenter CE if not already triggered from pre-timestep.
        for (size_t j = 1; j < Nactive; j++){
            if (_switch_peri(r, j)){
                ri_trace->current_C = 1;
                new_close_encounter = 1;
                if (ri_trace->peri_mode == REB_TRACE_PERI_FULL_BS || ri_trace->peri_mode == REB_TRACE_PERI_FULL_IAS15){
                    // Everything will be integrated with BS/IAS15. No need to check any further.
                    return new_close_encounter;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    ri_trace->tponly_encounter = 0;
                    break; // No need to check other particles
                }
            }
        }
    }
    if (ri_trace->current_C){
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        ri_trace->encounter_N = N;
        for (size_t i = 0; i < N; i++){
            ri_trace->encounter_map[i] = 1; // trigger encounter
        }
    }


    // Body-body
    // there cannot be TP-TP CEs
    for (size_t i = 0; i < Nactive; i++){ // Do not check for central body anymore
        for (size_t j = i + 1; j < N; j++){
            if (_switch(r, i, j)){
                if (ri_trace->current_Ks[i*N+j] == 0){
                    new_close_encounter = 1;
                }
                ri_trace->current_Ks[i*N+j] = 1;
                if (ri_trace->encounter_map[i] == 0){
                    ri_trace->encounter_map[i] = 1; // trigger encounter
                    ri_trace->encounter_N++;
                }
                if (ri_trace->encounter_map[j] == 0){
                    ri_trace->encounter_map[j] = 1; // trigger encounter
                    ri_trace->encounter_N++;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    ri_trace->tponly_encounter = 0;
                }
            }
        }
    }

    return new_close_encounter;
}

static void reb_integrator_trace_step_try(struct reb_simulation* const r){
    if (r->ri_trace.current_C == 0 || r->ri_trace.peri_mode == REB_TRACE_PERI_PARTIAL_BS){
        reb_integrator_trace_interaction_step(r, r->dt/2.);
        reb_integrator_trace_jump_step(r, r->dt/2.);
        reb_integrator_trace_kepler_step(r, r->dt);
        reb_integrator_trace_com_step(r,r->dt);
        reb_integrator_trace_jump_step(r, r->dt/2.);
        reb_integrator_trace_interaction_step(r, r->dt/2.);
    }else{
        // Pericenter approach with one of the FULL prescriptions
        double t_needed = r->t + r->dt;
        const double old_dt = r->dt;
        const double old_t = r->t;
        r->gravity = REB_GRAVITY_BASIC;
        r->ri_trace.mode = REB_TRACE_MODE_FULL;
        reb_integrator_trace_dh_to_inertial(r);
        switch (r->ri_trace.peri_mode){
            case REB_TRACE_PERI_FULL_IAS15:
                // Run default IAS15 integration
                reb_integrator_ias15_reset(r);
                while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 && r->status<=0){
                    reb_integrator_ias15_step(r);
                    if (r->t+r->dt >  t_needed){
                        r->dt = t_needed-r->t;
                    }
                    reb_collision_search(r);
                    if (r->N_collisions) r->ri_trace.force_accept = 1;
                }
                // Resetting IAS15 here reduces binary file size.
                reb_integrator_ias15_reset(r);
                break;
            case REB_TRACE_PERI_FULL_BS:
                {
                    // Run default BS integration
                    // TODO: Syntax should be similar to IAS
                    struct reb_ode* nbody_ode = NULL;

                    double* y;
                    while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 && r->status<=0){
                        if (!nbody_ode || nbody_ode->length != 6*r->N){
                            if (nbody_ode){
                                reb_ode_free(nbody_ode);
                            }
                            nbody_ode = reb_ode_create(r, 6*r->N);
                            nbody_ode->derivatives = reb_integrator_bs_nbody_derivatives;
                            nbody_ode->needs_nbody = 0;
                            y = nbody_ode->y;
                            reb_integrator_bs_reset(r);
                        }

                        for (size_t i=0; i<r->N; i++){
                            const struct reb_particle p = r->particles[i];
                            y[i*6+0] = p.x;
                            y[i*6+1] = p.y;
                            y[i*6+2] = p.z;
                            y[i*6+3] = p.vx;
                            y[i*6+4] = p.vy;
                            y[i*6+5] = p.vz;
                        }

                        int success = reb_integrator_bs_step_odes(r, r->dt);
                        if (success){
                            r->t += r->dt;
                        }
                        r->dt = r->ri_bs.dt_proposed;
                        if (r->t+r->dt >  t_needed){
                            r->dt = t_needed-r->t;
                        }

                        reb_integrator_bs_update_particles(r, nbody_ode->y);

                        if (success){
                            // Only do a collision search for accepted steps.
                            reb_collision_search(r);
                            if (r->N_collisions) r->ri_trace.force_accept = 1;
                        }
                    }
                    reb_ode_free(nbody_ode);
                    // Resetting BS here reduces binary file size
                    reb_integrator_bs_reset(r);
                }
                break;
            default:
                reb_simulation_error(r,"Unsupported peri_mode encountered\n");
                break;
        }
        r->t = old_t; // final time will be set later
        r->dt = old_dt;
        reb_integrator_trace_inertial_to_dh(r);
    }
}

void reb_integrator_trace_did_add_particle(struct reb_simulation* r){
    // TRACE can add particles mid-timestep now
    struct reb_integrator_trace* ri_trace = &(r->ri_trace);
    if (ri_trace->mode==REB_TRACE_MODE_KEPLER){
        const size_t old_N = r->N-1;
        if (ri_trace->N_allocated < r->N){
            ri_trace->current_Ks    = realloc(ri_trace->current_Ks, sizeof(int)*r->N*r->N);
            ri_trace->particles_backup = realloc(ri_trace->particles_backup, sizeof(struct reb_particle)*r->N);
            ri_trace->particles_backup_kepler = realloc(ri_trace->particles_backup_kepler, sizeof(struct reb_particle)*r->N);
            ri_trace->current_Ks    = realloc(ri_trace->current_Ks, sizeof(int)*r->N*r->N);
            ri_trace->encounter_map = realloc(ri_trace->encounter_map, sizeof(size_t)*r->N);
            r->map = ri_trace->encounter_map;
            ri_trace->encounter_map_backup = realloc(ri_trace->encounter_map_backup, sizeof(size_t)*r->N);
            ri_trace->N_allocated   = r->N;
        }

        // First reshuffle existing Ks
        size_t i = old_N;
        while (i --> 0){
            size_t j = old_N;
            while (j --> 0){
                ri_trace->current_Ks[i*old_N+j+i] = ri_trace->current_Ks[i*old_N+j];
            }
        }

        // add in new particle, we want it to interact with all currently interacting particles
        // exclude star
        for (size_t i = 1; i < ri_trace->encounter_N; i++){
            ri_trace->current_Ks[ri_trace->encounter_map[i]*r->N+old_N] = 1;
        }

        ri_trace->encounter_map[ri_trace->encounter_N] = old_N;
        ri_trace->encounter_N++;
        r->N_map++;

        if (r->N_active==SIZE_MAX){ 
            // If global N_active is not set, then all particles are active, so the new one as well.
            // Otherwise, assume we're adding non active particle. 
            ri_trace->encounter_N_active++;
        }

    }
}

void reb_integrator_trace_will_remove_particle(struct reb_simulation* r, size_t index){
    struct reb_integrator_trace* ri_trace = &(r->ri_trace);
    reb_integrator_bs_reset(r);
    if (ri_trace->mode==REB_TRACE_MODE_KEPLER){
        // Only removed mid-timestep if collision - BS Step!
        int after_to_be_removed_particle = 0;
        size_t encounter_index = SIZE_MAX;
        for (size_t i=0;i<ri_trace->encounter_N;i++){
            if (after_to_be_removed_particle == 1){
                ri_trace->encounter_map[i-1] = ri_trace->encounter_map[i] - 1;
            }
            if (ri_trace->encounter_map[i]==index){
                encounter_index = i;
                after_to_be_removed_particle = 1;
            }
        }
        if (encounter_index == SIZE_MAX){
            reb_simulation_error(r,"Cannot find particle in encounter map. Did not remove particle.");
            return;
        }

        // reshuffle current_Ks
        size_t counter = 0;
        const size_t new_N = r->N-1;
        for (size_t i = 0; i < new_N; i++){
            if (i == index) counter += r->N;
            for (size_t j = 0; j < new_N; j++){
                if (j == index) counter++;
                ri_trace->current_Ks[i*new_N+j] = ri_trace->current_Ks[i*new_N+j+counter];
            }
        }
        if (encounter_index<ri_trace->encounter_N_active){
            ri_trace->encounter_N_active--;
        }
        ri_trace->encounter_N--;
        r->N_map--;
    }
}


void reb_integrator_trace_step(struct reb_simulation* r){
    // Do memory management and consistency checks
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const size_t N = r->N;

    if (r->N_var){
        reb_simulation_warning(r,"TRACE does not work with variational equations.");
    }

    if (ri_trace->N_allocated<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_trace->particles_backup       = realloc(ri_trace->particles_backup,sizeof(struct reb_particle)*N);
        ri_trace->particles_backup_kepler   = realloc(ri_trace->particles_backup_kepler,sizeof(struct reb_particle)*N);
        ri_trace->current_Ks             = realloc(ri_trace->current_Ks,sizeof(int)*N*N);
        ri_trace->encounter_map          = realloc(ri_trace->encounter_map,sizeof(size_t)*N);
        ri_trace->encounter_map_backup   = realloc(ri_trace->encounter_map_backup,sizeof(size_t)*N);
        ri_trace->N_allocated = N;
    }

    // Calculate collisions only with DIRECT or LINE method
    if (r->collision != REB_COLLISION_NONE && (r->collision != REB_COLLISION_DIRECT && r->collision != REB_COLLISION_LINE)){
        reb_simulation_warning(r,"TRACE only works with a direct or line collision search.");
    }
    r->N_targets = SIZE_MAX; // Search for collisions between all particles in encounter step or full steps.

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_CUSTOM){
        reb_simulation_warning(r,"TRACE has its own gravity routine. Gravity routine set by the user will be ignored.");
    }

    // Not sure why this was needed. HR 4 April 2026
    // reb_integrator_trace_update_acceleration(r);

    reb_integrator_trace_inertial_to_dh(r);

    // Create copy of all particle to allow for the step to be rejected.
    memcpy(ri_trace->particles_backup, r->particles, N*sizeof(struct reb_particle));

    // This will be set to 1 if a collision occurred.
    ri_trace->force_accept = 0;

    // Check if there are any close encounters
    reb_integrator_trace_pre_ts_check(r);

    // Attempt one step. 
    reb_integrator_trace_step_try(r);

    // We always accept the step if a collision occurred as it is impossible to undo the collision.
    if (!ri_trace->force_accept){
        // We check again for close encounters to ensure time reversibility. 
        if (reb_integrator_trace_post_ts_check(r)){
            // New encounters were found. Will reject the step.
            // Revert particles to the beginning of the step.
            memcpy(r->particles, ri_trace->particles_backup, N*sizeof(struct reb_particle));

            // Do step again
            reb_integrator_trace_step_try(r);
        }
    }
    reb_integrator_trace_dh_to_inertial(r);

    r->t+=r->dt;
    r->dt_last_done = r->dt;
    r->N_targets = 1; // Only serch for collisions with star after complete timestep.
}

void reb_integrator_trace_reset(struct reb_simulation* r){
    r->ri_trace.encounter_N = 0;
    r->ri_trace.encounter_N_active = 0;
    r->ri_trace.r_crit_hill = 3;
    r->ri_trace.peri_crit_eta = 1.0;
    r->ri_trace.force_accept = 0;

    // Internal arrays (only used within one timestep)
    free(r->ri_trace.particles_backup);
    r->ri_trace.particles_backup = NULL;
    free(r->ri_trace.particles_backup_kepler);
    r->ri_trace.particles_backup_kepler = NULL;
    free(r->ri_trace.particles_backup_additional_forces);
    r->ri_trace.particles_backup_additional_forces = NULL;

    free(r->ri_trace.encounter_map);
    r->ri_trace.encounter_map = NULL;
    free(r->ri_trace.encounter_map_backup);
    r->ri_trace.encounter_map_backup = NULL;

    r->ri_trace.current_C = 0;
    free(r->ri_trace.current_Ks);
    r->ri_trace.current_Ks = NULL;

    r->ri_trace.S = NULL;
    r->ri_trace.S_peri = NULL;

    r->ri_trace.peri_mode = REB_TRACE_PERI_FULL_BS;

    r->ri_trace.N_allocated = 0;
    r->ri_trace.N_allocated_additional_forces = 0;
}
