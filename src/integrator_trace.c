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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_trace.h"
#include "integrator_whfast.h"
#include "integrator_bs.h"
#include "integrator_ias15.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

int reb_integrator_trace_switch_default(struct reb_simulation* const r, const unsigned int i, const unsigned int j){
    // for this test hard code no binary close encounter
    if (r->ri_trace.coordinates == REB_TRACE_COORDINATES_DHC && (j == 7 || i == 7)) return 0;

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

int reb_integrator_trace_switch_peri_default(struct reb_simulation* const r, const unsigned int j){
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

int reb_integrator_trace_switch_peri_none(struct reb_simulation* const r, const unsigned int j){
    // No pericenter flags
    return 0;
}

void reb_integrator_trace_inertial_to_dh(struct reb_simulation* r){
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
    // Particle 0 is also changed to allow for easy collision detection
    for (int i=N-1;i>=0;i--){
        particles[i].x -= particles[0].x;
        particles[i].y -= particles[0].y;
        particles[i].z -= particles[0].z;
        particles[i].vx -= com_vel.x;
        particles[i].vy -= com_vel.y;
        particles[i].vz -= com_vel.z;
    }
    r->ri_trace.com_pos = com_pos;
    r->ri_trace.com_vel = com_vel;
}

void reb_integrator_trace_inertial_to_wb(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    const int N = r->N;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;

    const int idxA = 0; // star A assumed to be first particle
    const int idxB = N_active-1; // star B assumed to be last active particle
    const int first_planet = 1; // planets index 1->N-2
    const int last_planet = N_active-2;

    double mplanets = 0.0; // total planet mass
    struct reb_vec3d sum_mx = {0.0, 0.0, 0.0};   // sum(m_i x_i)
    struct reb_vec3d sum_mvx = {0.0, 0.0, 0.0};  // sum(m_i v_i)

    for (int i = first_planet; i <= last_planet; ++i){
        double m = particles[i].m;
        mplanets += m;
        sum_mx.x  += m * particles[i].x;
        sum_mx.y  += m * particles[i].y;
        sum_mx.z  += m * particles[i].z;
        sum_mvx.x += m * particles[i].vx;
        sum_mvx.y += m * particles[i].vy;
        sum_mvx.z += m * particles[i].vz;
    }

    double mA = particles[idxA].m;
    double mB = particles[idxB].m;
    double mtot = mA + mB + mplanets;

    // ------------------ POSITIONS ------------------
    // Save B original position
    struct reb_vec3d xB_orig = {
        .x=particles[idxB].x,
        .y=particles[idxB].y,
        .z=particles[idxB].z
    };

    // X_i = x_i - x_A
    // Particle 0 is also changed to allow for easy collision detection NOT YET BUT PUT THIS BACK IN!!!!!
    for (int i=N-1; i>0; i--){
        if (i == idxB) continue; // skip binary companion
        particles[i].x -= particles[idxA].x;
        particles[i].y -= particles[idxA].y;
        particles[i].z -= particles[idxA].z;
    }

    // X_B = x_B - (mA xA + sum_mx)/(mA+mplanets)
    particles[idxB].x -= (mA*particles[idxA].x + sum_mx.x)/(mA+mplanets);
    particles[idxB].y -= (mA*particles[idxA].y + sum_mx.y)/(mA+mplanets);
    particles[idxB].z -= (mA*particles[idxA].z + sum_mx.z)/(mA+mplanets);

    // X_A = (mA xA + mB xB + sum_mx) / mtot
    particles[idxA].x = (mA * particles[idxA].x + mB * xB_orig.x + sum_mx.x) / mtot;
    particles[idxA].y = (mA * particles[idxA].y + mB * xB_orig.y + sum_mx.y) / mtot;
    particles[idxA].z = (mA * particles[idxA].z + mB * xB_orig.z + sum_mx.z) / mtot;

    // ------------------ VELOCITIES ------------------
    // Save original inertial velocities before any writes
    struct reb_vec3d vA_orig = { particles[idxA].vx, particles[idxA].vy, particles[idxA].vz };
    struct reb_vec3d vB_orig = { particles[idxB].vx, particles[idxB].vy, particles[idxB].vz };

    // Planets: V_i = v_i - (mA vA + sum_mvx)/(mA+mplanets)
    // Particle 0 is also changed to allow for easy collision detection NOT YET BUT PUT THIS BACK IN!!!!!
    for (int i=N-1; i>0; i--) {
        if (i == idxB) continue; // skip binary companion
        particles[i].vx -= (mA*vA_orig.x + sum_mvx.x)/(mA+mplanets);
        particles[i].vy -= (mA*vA_orig.y + sum_mvx.y)/(mA+mplanets);
        particles[i].vz -= (mA*vA_orig.z + sum_mvx.z)/(mA+mplanets);
    }

    // Total momentum (inertial)
    struct reb_vec3d Ptot = {
        .x = mA*vA_orig.x + mB*vB_orig.x + sum_mvx.x,
        .y = mA*vA_orig.y + mB*vB_orig.y + sum_mvx.y,
        .z = mA*vA_orig.z + mB*vB_orig.z + sum_mvx.z
    };

    // A DOF is the COM: V_A = Ptot / mtot
    particles[idxA].vx = Ptot.x / mtot;
    particles[idxA].vy = Ptot.y / mtot;
    particles[idxA].vz = Ptot.z / mtot;

    // B DOF: V_B = v_B - Ptot/mtot  (i.e. v_B - V_A)
    particles[idxB].vx -= Ptot.x / mtot;
    particles[idxB].vy -= Ptot.y / mtot;
    particles[idxB].vz -= Ptot.z / mtot;

    // Save COM info (inertial COM before transform)
    r->ri_trace.com_pos.x = particles[idxA].x;
    r->ri_trace.com_pos.y = particles[idxA].y;
    r->ri_trace.com_pos.z = particles[idxA].z;
    r->ri_trace.com_vel.x = Ptot.x / mtot;
    r->ri_trace.com_vel.y = Ptot.y / mtot;
    r->ri_trace.com_vel.z = Ptot.z / mtot;
}

// CHATGPT function, seems to work correctly but need to carefully verify
void reb_integrator_trace_wb_to_inertial(struct reb_simulation* r) {
    struct reb_particle* particles = r->particles;
    const int N = r->N;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;

    const int idxA = 0;        // "A DOF": stores X_A, V_A
    const int idxB = N_active - 1;    // "B DOF": stores X_B, V_B
    const int first_planet = 1;
    const int last_planet  = N_active - 2;

    // ------------------ sums from WB planet coordinates ------------------
    double M = 0.0; // Σ m_i
    struct reb_vec3d sum_mX = {0.0, 0.0, 0.0}; // Σ m_i X_i
    struct reb_vec3d sum_mV = {0.0, 0.0, 0.0}; // Σ m_i V_i

    for (int i = first_planet; i <= last_planet; ++i) {
        const double m = particles[i].m;
        M += m;
        sum_mX.x += m * particles[i].x;
        sum_mX.y += m * particles[i].y;
        sum_mX.z += m * particles[i].z;

        sum_mV.x += m * particles[i].vx;
        sum_mV.y += m * particles[i].vy;
        sum_mV.z += m * particles[i].vz;
    }

    const double mA = particles[idxA].m;
    const double mB = particles[idxB].m;
    const double mtot   = mA + mB + M;
    const double denomA = mA + M;  // (mA + Σ m_i)

    // Cache WB star DOFs (avoid clobbering)
    const struct reb_vec3d X_A = { particles[idxA].x,  particles[idxA].y,  particles[idxA].z };
    const struct reb_vec3d X_B = { particles[idxB].x,  particles[idxB].y,  particles[idxB].z };
    const struct reb_vec3d V_A = { particles[idxA].vx, particles[idxA].vy, particles[idxA].vz }; // COM velocity
    const struct reb_vec3d V_B = { particles[idxB].vx, particles[idxB].vy, particles[idxB].vz };

    // ------------------ POSITIONS inverse ------------------
    // Using:
    //   X_i = x_i - x_A                      => x_i = X_i + x_A
    //   X_B = x_B - (mA x_A + Σ m_i x_i)/(mA+M)
    //   X_A = (mA x_A + mB x_B + Σ m_i x_i)/m_tot
    //
    // Let Σ m_i x_i = Σ m_i (X_i + x_A) = Σ m_i X_i + M x_A = sum_mX + M x_A
    // Then (mA x_A + Σ m_i x_i) = (mA+M) x_A + sum_mX = denomA * x_A + sum_mX
    // So:
    //   x_B = X_B + x_A + sum_mX/denomA
    // Plug into X_A:
    //   m_tot X_A = mA x_A + mB x_B + Σ m_i x_i
    //            = m_tot x_A + mB X_B + mB (sum_mX/denomA) + sum_mX
    // => x_A = X_A - [ mB X_B + mB (sum_mX/denomA) + sum_mX ] / m_tot

    struct reb_vec3d sum_mX_over_d = {0.0, 0.0, 0.0};
    sum_mX_over_d.x = sum_mX.x / denomA;
    sum_mX_over_d.y = sum_mX.y / denomA;
    sum_mX_over_d.z = sum_mX.z / denomA;

    struct reb_vec3d xA = X_A;
    xA.x = X_A.x - (mB * X_B.x + mB * sum_mX_over_d.x + sum_mX.x) / mtot;
    xA.y = X_A.y - (mB * X_B.y + mB * sum_mX_over_d.y + sum_mX.y) / mtot;
    xA.z = X_A.z - (mB * X_B.z + mB * sum_mX_over_d.z + sum_mX.z) / mtot;

    // Planets: x_i = X_i + x_A
    for (int i=1;i<N;i++) {
        if (i == idxB) continue; // skip binary companion
        particles[i].x += xA.x; // particles[i].x currently X_i
        particles[i].y += xA.y;
        particles[i].z += xA.z;
    }

    // Star B: x_B = X_B + x_A + sum_mX/denomA
    const struct reb_vec3d xB = {
        X_B.x + xA.x + sum_mX_over_d.x,
        X_B.y + xA.y + sum_mX_over_d.y,
        X_B.z + xA.z + sum_mX_over_d.z
    };

    // Write inertial star positions
    particles[idxA].x = xA.x; particles[idxA].y = xA.y; particles[idxA].z = xA.z;
    particles[idxB].x = xB.x; particles[idxB].y = xB.y; particles[idxB].z = xB.z;

    // ------------------ VELOCITIES inverse ------------------
    // Forward (with your corrected V_A):
    //   V_A = P_tot / m_tot  (COM velocity)
    //   V_B = v_B - V_A
    //   V_i = v_i - W, where W = (mA v_A + Σ m_i v_i)/(mA+M)
    //
    // Invert:
    //   v_B = V_B + V_A
    //   W   = V_A - (mB/(mA+M)) V_B      (from COM momentum algebra)
    //   v_i = V_i + W = V_i + V_A - (mB/(mA+M)) V_B
    //   v_A = V_A - (mB/(mA+M)) V_B - (1/mA) Σ m_i V_i

    const double coefB = (denomA != 0.0) ? (mB / denomA) : 0.0;

    // Common additive term for all planets: add = V_A - coefB * V_B
    const struct reb_vec3d add_planet = {
        V_A.x - coefB * V_B.x,
        V_A.y - coefB * V_B.y,
        V_A.z - coefB * V_B.z
    };

    // v_A
    struct reb_vec3d vA = {
        add_planet.x - ((mA != 0.0) ? (sum_mV.x / mA) : 0.0),
        add_planet.y - ((mA != 0.0) ? (sum_mV.y / mA) : 0.0),
        add_planet.z - ((mA != 0.0) ? (sum_mV.z / mA) : 0.0)
    };

    // v_B
    struct reb_vec3d vB = {
        V_B.x + V_A.x,
        V_B.y + V_A.y,
        V_B.z + V_A.z
    };

    // Planets: v_i = V_i + add_planet
    for (int i = first_planet; i <= N; ++i) {
        if (i == idxB) continue;
        particles[i].vx += add_planet.x; // particles[i].vx currently V_i
        particles[i].vy += add_planet.y;
        particles[i].vz += add_planet.z;
    }

    // Write inertial star velocities
    particles[idxA].vx = vA.x; particles[idxA].vy = vA.y; particles[idxA].vz = vA.z;
    particles[idxB].vx = vB.x; particles[idxB].vy = vB.y; particles[idxB].vz = vB.z;
}


void reb_integrator_trace_dh_to_inertial(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle temp = {0};
    const int N = r->N;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;
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
    // Should not matter unless collisions occurred.
    particles[0].x = r->ri_trace.com_pos.x - temp.x;
    particles[0].y = r->ri_trace.com_pos.y - temp.y;
    particles[0].z = r->ri_trace.com_pos.z - temp.z;

    for (int i=1;i<N;i++){
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

void reb_integrator_trace_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    r->ri_trace.mode = REB_TRACE_MODE_INTERACTION;
    reb_simulation_update_acceleration(r);
    for (int i=1;i<N;i++){
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

    const int N_active = r->N_active==-1?r->N:r->N_active;
    const int N_all = r->N;

    // If TP type 1, use r->N. Else, use N_active.
    int N = r->testparticle_type==0 ? N_active: r->N;

    // If wb coordinates, don't include star B
    const int idxB = N_active-1;

    double px=0., py=0., pz=0.;
    for (int i=1;i<N;i++){
        if (r->ri_trace.coordinates==REB_TRACE_COORDINATES_WB && i==idxB) continue; // skip star B
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m;
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px *= dt/r->particles[0].m;
    py *= dt/r->particles[0].m;
    pz *= dt/r->particles[0].m;

    for (int i=1;i<N_all;i++){
        if (r->ri_trace.coordinates==REB_TRACE_COORDINATES_WB && i==idxB) continue; // skip star B
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

/*
void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    //struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    double total_mass = r->particles[0].m;
    for (int i=1;i<N;i++){
        // if WB, we need to change the central mass
        if (i==N-1 && r->ri_trace.coordinates==REB_TRACE_COORDINATES_WB){
            reb_whfast_kepler_solver(r,r->particles,r->G*total_mass,i,dt);
            break;
        }
        else reb_whfast_kepler_solver(r,r->particles,r->G*r->particles[0].m,i,dt);
        total_mass += r->particles[i].m;
    }
}
    */

void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const int idxB = N_active-1;
    struct reb_particle* p = r->particles;

    double mA = p[0].m;
    double mB = p[idxB].m;
    double Mpl = 0.0;
    for (int i=1; i<N; ++i){
        if (r->ri_trace.coordinates==REB_TRACE_COORDINATES_WB && i==idxB){
            // Kepler parameter for relative coordinate X_B is G * m_tot
            const double GM = r->G * (mA + Mpl + mB);

            // I believe we need to convert from V_B equiv P_B / m_B to real velocities here, just for the Kepler solver.
            // Works! note to self get this sanity checked...

            const double scale = ((mA + mB + Mpl) / (mA + Mpl));
            p[i].vx *= scale;  p[i].vy *= scale;  p[i].vz *= scale;
            reb_whfast_kepler_solver(r, p, GM, i, dt);
            p[i].vx /= scale;  p[i].vy /= scale;  p[i].vz /= scale;
        }else{
            // normal Keplerian orbit
            reb_whfast_kepler_solver(r, p, r->G*p[0].m, i, dt);
            Mpl += p[i].m;
        }
    }
}


void reb_integrator_trace_update_particles(struct reb_simulation* r, const double* y){
    int N = r->ri_trace.encounter_N;
    int* map = r->ri_trace.encounter_map;

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

void reb_integrator_trace_nbody_derivatives(struct reb_ode* ode, double* const yDot, const double* const y, double const t){
    struct reb_simulation* const r = ode->r;
    // TRACE always needs this to ensure the right Hamiltonian is evolved
    reb_integrator_trace_update_particles(r, y);
    reb_simulation_update_acceleration(r);

    double px=0., py=0., pz=0.;
    int* map = r->ri_trace.encounter_map;
    int N = r->ri_trace.encounter_N;

    if (map==NULL){
        reb_simulation_error(r, "Cannot access TRACE map from BS.");
        return;
    }

    // Kepler Step
    // This is only for pericenter approach
    if (r->ri_trace.current_C){
        for (int i=1;i<r->N;i++){ // all particles
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

    for (int i=1; i<N; i++){
        int mi = map[i];
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

    int i_enc = 0;
    ri_trace->encounter_N_active = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(ri_trace->encounter_map[i]){
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = ri_trace->particles_backup_kepler[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            ri_trace->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
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

    if (ri_trace->peri_mode == REB_TRACE_PERI_PARTIAL_BS || !ri_trace->current_C){
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

            // Initialize ODE
            for (unsigned int i=0; i<ri_trace->encounter_N; i++){
                const int mi = ri_trace->encounter_map[i];
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
                if (r->collisions_N) r->ri_trace.force_accept = 1;
            }

            if (nbody_ode->length != ri_trace->encounter_N*3*2){
                // Just re-create the ODE
                reb_ode_free(nbody_ode);
                nbody_ode = reb_ode_create(r, ri_trace->encounter_N*3*2);
                nbody_ode->derivatives = reb_integrator_trace_nbody_derivatives;
                nbody_ode->needs_nbody = 0;
                r->ri_bs.first_or_last_step = 1;
            }

            star.vx = r->particles[0].vx; // keep track of changed star velocity for later collisions
            star.vy = r->particles[0].vy;
            star.vz = r->particles[0].vz;

            if (r->particles[0].x !=0 || r->particles[0].y !=0 || r->particles[0].z !=0){
                // Collision with star occurred
                // Shift all particles back to heliocentric coordinates
                // Ignore stars velocity:
                //   - will not be used after this
                //   - com velocity is unchanged. this velocity will be used
                //     to reconstruct star's velocity later.
                for (int i=r->N-1; i>=0; i--){
                    r->particles[i].x -= r->particles[0].x;
                    r->particles[i].y -= r->particles[0].y;
                    r->particles[i].z -= r->particles[0].z;
                }
            }
        }

        // if only test particles encountered massive bodies, reset the
        // massive body coordinates to their post Kepler step state
        if(ri_trace->tponly_encounter){
            for (unsigned int i=1; i < ri_trace->encounter_N_active; i++){
                unsigned int mi = ri_trace->encounter_map[i];
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
}

void reb_integrator_trace_kepler_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    memcpy(ri_trace->particles_backup_kepler,r->particles,r->N*sizeof(struct reb_particle));
    reb_integrator_trace_whfast_step(r, _dt);
    reb_integrator_trace_bs_step(r, _dt);
}


void reb_integrator_trace_pre_ts_check(struct reb_simulation* const r){
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;
    int idxB = Nactive - 1; // only used for WB
    int (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = ri_trace->S ? ri_trace->S : reb_integrator_trace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = ri_trace->S_peri ? ri_trace->S_peri : reb_integrator_trace_switch_peri_default;

    // Clear encounter map
    for (unsigned int i=1; i<r->N; i++){
        ri_trace->encounter_map[i] = 0;
    }
    ri_trace->encounter_map[0] = 1;
    ri_trace->encounter_N = 1;

    // Reset encounter triggers.
    ri_trace->current_C = 0;

    for (int i = 0; i < r->N; i++){
        for (unsigned int j = i + 1; j < r->N; j++){
            ri_trace->current_Ks[i*r->N+j] = 0;
        }
    }

    if (r->testparticle_type == 1){
        ri_trace->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        ri_trace->tponly_encounter = 1;
    }

    // Don't think this needs to change for WB
    // Check for pericenter CE
    for (int j = 1; j < Nactive; j++){
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
        for (int i = 1; i < N; i++){
            ri_trace->encounter_map[i] = 1; //  trigger encounter
        }

    }

    // Body-body
    // there cannot be TP-TP CEs
    for (int i = 0; i < Nactive; i++){ // Check central body, for collisions
        for (int j = i + 1; j < N; j++){
            if ((i == idxB || j == idxB) && ri_trace->coordinates == REB_TRACE_COORDINATES_WB) continue; // in WB coordinates star B does not have close encounters, for now
            if (_switch(r, i, j)){
                ri_trace->current_Ks[i*r->N+j] = 1;
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
    memcpy(ri_trace->encounter_map_backup, ri_trace->encounter_map, N*sizeof(int));
}

double reb_integrator_trace_post_ts_check(struct reb_simulation* const r){
    // This function returns 1 if any new encounters occurred.
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;
    const int Nactive = r->N_active==-1?r->N:r->N_active;
    int idxB = Nactive - 1; // only used for WB
    int (*_switch) (struct reb_simulation* const r, const unsigned int i, const unsigned int j) = ri_trace->S ? ri_trace->S : reb_integrator_trace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const unsigned int j) = ri_trace->S_peri ? ri_trace->S_peri : reb_integrator_trace_switch_peri_default;
    int new_close_encounter = 0; // New CEs

    // Clear encounter maps
    /*
       for (unsigned int i=1; i<r->N; i++){
       ri_trace->encounter_map[i] = 0;
       }
       ri_trace->encounter_map[0] = 1;
       ri_trace->encounter_N = 1;
     */

    // Set this from pre-ts encounter map. I don't think we need to reset encounter_N here.
    memcpy(ri_trace->encounter_map, ri_trace->encounter_map_backup, N*sizeof(int));

    if (!ri_trace->current_C){
        // Check for pericenter CE if not already triggered from pre-timestep.
        for (int j = 1; j < Nactive; j++){
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
        for (int i = 0; i < N; i++){
            ri_trace->encounter_map[i] = 1; // trigger encounter
        }
    }


    // Body-body
    // there cannot be TP-TP CEs
    for (int i = 0; i < Nactive; i++){ // Do not check for central body anymore
        for (int j = i + 1; j < N; j++){
            if ((i == idxB || j == idxB) && ri_trace->coordinates == REB_TRACE_COORDINATES_WB) continue; // in WB coordinates star B does not have close encounters, for now
            if (_switch(r, i, j)){
                if (ri_trace->current_Ks[i*r->N+j] == 0){
                    new_close_encounter = 1;
                }
                ri_trace->current_Ks[i*r->N+j] = 1;
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

// TODO: Should be reused from BS
static void nbody_derivatives(struct reb_ode* ode, double* const yDot, const double* const y, double const t){
    struct reb_simulation* const r = ode->r;
    reb_integrator_bs_update_particles(r, y);
    reb_simulation_update_acceleration(r);

    for (unsigned int i=0; i<r->N; i++){
        const struct reb_particle p = r->particles[i];
        yDot[i*6+0] = p.vx;
        yDot[i*6+1] = p.vy;
        yDot[i*6+2] = p.vz;
        yDot[i*6+3] = p.ax;
        yDot[i*6+4] = p.ay;
        yDot[i*6+5] = p.az;
    }
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
        r->ri_trace.mode = REB_TRACE_MODE_FULL; // for collision search
        if (r->ri_trace.coordinates == REB_TRACE_COORDINATES_WB) reb_integrator_trace_wb_to_inertial(r);
        else reb_integrator_trace_dh_to_inertial(r);
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
                    if (r->collisions_N) r->ri_trace.force_accept = 1;
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
                            nbody_ode->derivatives = nbody_derivatives;
                            nbody_ode->needs_nbody = 0;
                            y = nbody_ode->y;
                            reb_integrator_bs_reset(r);
                        }

                        for (unsigned int i=0; i<r->N; i++){
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
                            if (r->collisions_N) r->ri_trace.force_accept = 1;
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
        if (r->ri_trace.coordinates == REB_TRACE_COORDINATES_WB) r->gravity = REB_GRAVITY_WB;
        else r->gravity = REB_GRAVITY_TRACE;
        r->t = old_t; // final time will be set later
        r->dt = old_dt;
        if (r->ri_trace.coordinates == REB_TRACE_COORDINATES_WB) reb_integrator_trace_inertial_to_wb(r);
        else reb_integrator_trace_inertial_to_dh(r);
    }
}

void reb_integrator_trace_step(struct reb_simulation* r){
    // Do memory management and consistency checks
    struct reb_integrator_trace* const ri_trace = &(r->ri_trace);
    const int N = r->N;

    if (r->N_var_config){
        reb_simulation_warning(r,"TRACE does not work with variational equations.");
    }

    if (ri_trace->N_allocated<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        ri_trace->particles_backup       = realloc(ri_trace->particles_backup,sizeof(struct reb_particle)*N);
        ri_trace->particles_backup_kepler   = realloc(ri_trace->particles_backup_kepler,sizeof(struct reb_particle)*N);
        ri_trace->current_Ks             = realloc(ri_trace->current_Ks,sizeof(int)*N*N);
        ri_trace->encounter_map          = realloc(ri_trace->encounter_map,sizeof(int)*N);
        ri_trace->encounter_map_backup   = realloc(ri_trace->encounter_map_backup,sizeof(int)*N);
        ri_trace->N_allocated = N;
    }

    // Calculate collisions only with DIRECT or LINE method
    if (r->collision != REB_COLLISION_NONE && (r->collision != REB_COLLISION_DIRECT && r->collision != REB_COLLISION_LINE)){
        reb_simulation_warning(r,"TRACE only works with a direct or line collision search.");
    }

    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_TRACE){
        reb_simulation_warning(r,"TRACE has its own gravity routine. Gravity routine set by the user will be ignored.");
    }

    if (ri_trace->coordinates == REB_TRACE_COORDINATES_WB) r->gravity = REB_GRAVITY_WB;
    else r->gravity = REB_GRAVITY_TRACE;

    ri_trace->mode = REB_TRACE_MODE_NONE; // Do not calculate gravity in-between timesteps. TRACE will call reb_update_acceleration itself.

    reb_simulation_update_acceleration(r);

    if (ri_trace->coordinates == REB_TRACE_COORDINATES_WB) reb_integrator_trace_inertial_to_wb(r);
    else reb_integrator_trace_inertial_to_dh(r);

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

    if (ri_trace->coordinates == REB_TRACE_COORDINATES_WB) reb_integrator_trace_wb_to_inertial(r);
    else reb_integrator_trace_dh_to_inertial(r);

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

void reb_integrator_trace_synchronize(struct reb_simulation* r){
}

void reb_integrator_trace_reset(struct reb_simulation* r){
    r->ri_trace.mode = REB_TRACE_MODE_NONE;
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
    r->ri_trace.coordinates = REB_TRACE_COORDINATES_DHC;

    r->ri_trace.N_allocated = 0;
    r->ri_trace.N_allocated_additional_forces = 0;
}
