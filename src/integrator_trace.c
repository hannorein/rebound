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
#include "binarydata.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

void reb_integrator_trace_step(struct reb_simulation* r, void* state);
void* reb_integrator_trace_create();
void reb_integrator_trace_free(void* p);
void reb_integrator_trace_did_add_particle(struct reb_simulation* r);
void reb_integrator_trace_will_remove_particle(struct reb_simulation* r, size_t index);
const struct reb_binarydata_field_descriptor reb_integrator_trace_field_descriptor_list[];

const struct reb_integrator reb_integrator_trace = {
    .documentation = 
    "TRACE is a hybrid (almost) time-reversible integrator, based on the algorithm described "
    "in [Hernandez & Dehnen (2023)]. It uses WHFast for long term integrations but "
    "switches to BS or IAS15 for all close encounters. TRACE is appropriate for systems "
    "with a dominant central mass if particles occasionally have close encounters. "
    "This TRACE implementation is described in [Lu, Hernandez & Rein (2024)]. " 
    "\n\n"
    "[Hernandez & Dehnen (2023)]: https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.4639H/abstract\n" 
    "[Lu, Hernandez & Rein (2024)]: https://ui.adsabs.harvard.edu/abs/2024MNRAS.533.3708L/abstract\n" 
    "[Pham, Rein, and Spiegel (2024)]: https://ui.adsabs.harvard.edu/abs/2024OJAp....7E...1P/abstract\n"
    ,
    .step = reb_integrator_trace_step,
    .create = reb_integrator_trace_create,
    .free = reb_integrator_trace_free,
    .will_remove_particle = reb_integrator_trace_will_remove_particle,
    .did_add_particle = reb_integrator_trace_did_add_particle,
    .field_descriptor_list = reb_integrator_trace_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_trace_field_descriptor_list[] = {
    { "The critical switchover radii of non-central particles are calculated based on a "
        "modified Hill radii criteria. This modified Hill radius for each particle is "
        "calculated and then multiplied by this parameter. The parameter is in units of "
        "the modified Hill radius. This value is used by the `default` switching function. "
        "The default value is 4.", 
        REB_DOUBLE,      "r_crit_hill",      offsetof(struct reb_integrator_trace_state, r_crit_hill), 0, 0, 0},
    { "Critical distance (units of hill radius) for close encounters with the binary in WB coordinates.", 
        REB_DOUBLE,      "r_crit_WB",      offsetof(struct reb_integrator_trace_state, r_crit_WB), 0, 0, 0},
    { "The criteria for a pericenter approach with the central body. This criteria is used "
        "in the `default` pericenter switching condition. It flags a particle as in a close "
        "pericenter approach if the ratio of the timestep to the condition described in "
        "[Pham, Rein, and Spiegel (2024)]. The default value is 1.",
        REB_DOUBLE,      "peri_crit_eta",    offsetof(struct reb_integrator_trace_state, peri_crit_eta), 0, 0, 0},
    { "This flag determines how TRACE integrates close approaches with the central star.", 
        REB_INT,         "peri_mode",        offsetof(struct reb_integrator_trace_state, peri_mode), 0, 0, REB_GENERATE_ENUM_DESCRIPTORS(REB_INTEGRATOR_TRACE_PERIMODE)},
    { "Coordinates system.", 
        REB_INT,         "coordinates",      offsetof(struct reb_integrator_trace_state, coordinates), 0, 0, REB_GENERATE_ENUM_DESCRIPTORS(REB_INTEGRATOR_TRACE_COORDINATES)},
    { "This is a function pointer to the switching function for close encounters between "
        "non-central bodies. If NULL (the default), the default switching function will be used."
        "The default switching function is a similar (but slightly modified) switching function "
        "used by MERCURY. It uses a modified Hill radius criteria, with heliocentric distance "
        "replacing the semimajor axis. "
        "\n\n"
        "The arguments `i` and `j` are the indices of the two particles considered. The return "
        "value is either 0 or 1. A return value of 1 means a close encounter has been flagged. "
        "If the return values of both this function and the central switching function are "
        "always 0, then the integrator effectively becomes the standard Wisdom-Holman integrator. ", 
        REB_FUNCTIONPOINTER,"S",            offsetof(struct reb_integrator_trace_state, S), 0, 0, 0},
    { "This is a function pointer to the switching function for close encounters involving "
        "the central body. If NULL (the default), the default switching function will be used. "
        "The default switching function checks if a body is close to its pericenter by "
        "considering a timescale derived from high-order derivatives of the particle's "
        "herliocentric position, inspired by [Pham, Rein, and Spiegel (2024)]. "
        "The argument `j` is the index of the non-central particle considered. "
        "The return value is either 0 or 1. A return value of 1 means a close encounter "
        "has been flagged. ",
        REB_FUNCTIONPOINTER,"S_peri",            offsetof(struct reb_integrator_trace_state, S), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* reb_integrator_trace_create(){
    struct reb_integrator_trace_state* trace = calloc(sizeof(struct reb_integrator_trace_state),1);
    trace->r_crit_hill = 3;
    trace->r_crit_WB = 0.1;
    trace->peri_crit_eta = 1.0;
    trace->S = NULL;
    trace->S_peri = NULL;
    trace->peri_mode = REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS;
    trace->coordinates = REB_INTEGRATOR_TRACE_COORDINATES_DEMOCRATICHELIOCENTRIC;
    return trace;
}

void reb_integrator_trace_free(void* state){
    struct reb_integrator_trace_state* trace = state;
    free(trace->particles_backup);
    free(trace->particles_backup_kepler);
    free(trace->particles_backup_additional_forces);
    free(trace->encounter_map);
    free(trace->encounter_map_backup);
    free(trace);
}

int reb_integrator_trace_switch_default(struct reb_simulation* const r, const size_t i, const size_t j){
    // Returns 1 for close encounter between i and j, 0 otherwise
    struct reb_integrator_trace_state* trace = r->integrator.state;
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

    double r_crit_hill2 = trace->r_crit_hill*trace->r_crit_hill;
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
    const struct reb_integrator_trace_state* const trace = r->integrator.state;
    double GM = r->G*r->particles[0].m; // Not sure if this is the right mass to use.

    double x = r->particles[j].x - r->particles[0].x;
    double y = r->particles[j].y - r->particles[0].y;
    double z = r->particles[j].z - r->particles[0].z;
    double d2 = x*x + y*y + z*z;
    double d = sqrt(d2);

    // first derivative
    double dx = r->particles[j].vx - r->particles[0].vx;
    double dy = r->particles[j].vy - r->particles[0].vy;
    double dz = r->particles[j].vz - r->particles[0].vz;

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
    double dt_prs2 = trace->peri_crit_eta * trace->peri_crit_eta * tau_prs2;

    if (r->dt * r->dt > dt_prs2){
        //printf("This happening??? %f\n", r->t);
        //exit(1);
        return 1;
    }else{
        // In WB coordinates, we also switch if we enter the binary's hill sphere x r_crit_WB

        if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY){
            // reb_integrator_trace_wb_to_inertial(r);
            const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary")); // star B assumed to be last active particle
            const int has_binary = (idxB != -1);
            if (has_binary){
                double bx = r->particles[idxB].x;
                double by = r->particles[idxB].y;
                double bz = r->particles[idxB].z;
                double br2 = bx*bx + by*by + bz*bz;
                double mr = r->particles[idxB].m/(3.*r->particles[0].m);

                double factor2 = trace->r_crit_WB*trace->r_crit_WB;
                const double rhillb6 = factor2*factor2*factor2*br2*br2*br2*mr*mr;

                double dxb = r->particles[j].x - r->particles[idxB].x;
                double dyb = r->particles[j].y - r->particles[idxB].y;
                double dzb = r->particles[j].z - r->particles[idxB].z;
                double d2b = dxb*dxb + dyb*dyb + dzb*dzb;

                //printf("this happening??? %f %d %d %f %f\n", r->t, j, idxB, sqrt(d2b), pow(rhillb6,1./6.));
                //    exit(1);

                if (d2b*d2b*d2b < rhillb6){
                    return 1;
                }
            }

        }
        return 0;
    }
}

int reb_integrator_trace_switch_WB_default(struct reb_simulation* const r, const size_t j){
    // reb_integrator_trace_wb_to_inertial(r);
    const struct reb_integrator_trace_state* const trace = r->integrator.state;
    const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary")); // star B assumed to be last active particle
    double bx = r->particles[idxB].x - r->particles[0].x;
    double by = r->particles[idxB].y - r->particles[0].y;
    double bz = r->particles[idxB].z - r->particles[0].z;
    double br2 = bx*bx + by*by + bz*bz;
    double mr = r->particles[idxB].m/(3.*r->particles[0].m);

    double factor2 = trace->r_crit_WB*trace->r_crit_WB;
    const double rhillb6 = factor2*factor2*factor2*br2*br2*br2*mr*mr;
    //const double rhillb62 = factor2*factor2*factor2*br2*br2*br2*mr*mr;

    double dxb = r->particles[j].x - r->particles[idxB].x;
    double dyb = r->particles[j].y - r->particles[idxB].y;
    double dzb = r->particles[j].z - r->particles[idxB].z;
    double d2b = dxb*dxb + dyb*dyb + dzb*dzb;

    if (d2b*d2b*d2b < rhillb6){
        return 1;
    }
    return 0;
}

int reb_integrator_trace_switch_peri_none(struct reb_simulation* const r, const size_t j){
    // No pericenter flags
    (void)r; // Not used
    (void)j; // Not used
    return 0;
}

void reb_integrator_trace_inertial_to_dh(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
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
    trace->com_pos = com_pos;
    trace->com_vel = com_vel;
}

void reb_integrator_trace_dh_to_inertial(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
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
    particles[0].x = trace->com_pos.x - temp.x;
    particles[0].y = trace->com_pos.y - temp.y;
    particles[0].z = trace->com_pos.z - temp.z;

    for (size_t i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += trace->com_vel.x;
        particles[i].vy += trace->com_vel.y;
        particles[i].vz += trace->com_vel.z;
    }
    particles[0].vx = trace->com_vel.x - temp.vx;
    particles[0].vy = trace->com_vel.y - temp.vy;
    particles[0].vz = trace->com_vel.z - temp.vz;
}

void reb_integrator_trace_inertial_to_wb(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const int N = r->N;
    const int N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?r->N:r->N_active;

    const int idxA = 0; // star A assumed to be first particle
    const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary"));

    if (idxB == -1){
        // no binary companion, just shift to DH coordinates
        reb_integrator_trace_inertial_to_dh(r);
        return;
    }

    double mplanets = 0.0; // total planet mass
    struct reb_vec3d sum_mx = {0.0, 0.0, 0.0};   // sum(m_i x_i)
    struct reb_vec3d sum_mvx = {0.0, 0.0, 0.0};  // sum(m_i v_i)

    for (int i = 1; i <= N_active-1; ++i){
        if (i == idxB) continue; // skip binary companion
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
    // X_A = COM position — save it, then zero particle A
    trace->com_pos.x = (mA * particles[idxA].x + mB * xB_orig.x + sum_mx.x) / mtot;
    trace->com_pos.y = (mA * particles[idxA].y + mB * xB_orig.y + sum_mx.y) / mtot;
    trace->com_pos.z = (mA * particles[idxA].z + mB * xB_orig.z + sum_mx.z) / mtot;
    // Particle 0 is also changed to allow for easy collision detection
    particles[idxA].x = 0.0;
    particles[idxA].y = 0.0;
    particles[idxA].z = 0.0;

    // ------------------ VELOCITIES ------------------
    // Save original inertial velocities before any writes
    struct reb_vec3d vA_orig = { particles[idxA].vx, particles[idxA].vy, particles[idxA].vz };
    struct reb_vec3d vB_orig = { particles[idxB].vx, particles[idxB].vy, particles[idxB].vz };

    // Planets: V_i = v_i - (mA vA + sum_mvx)/(mA+mplanets)
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
    trace->com_vel.x = Ptot.x / mtot;
    trace->com_vel.y = Ptot.y / mtot;
    trace->com_vel.z = Ptot.z / mtot;
    particles[idxA].vx = 0.0;
    particles[idxA].vy = 0.0;
    particles[idxA].vz = 0.0;

    // B DOF: V_B = v_B - Ptot/mtot  (i.e. v_B - V_A)
    particles[idxB].vx -= trace->com_vel.x;
    particles[idxB].vy -= trace->com_vel.y;
    particles[idxB].vz -= trace->com_vel.z;
}

void reb_integrator_trace_wb_to_inertial(struct reb_simulation* r) {
    struct reb_particle* particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const int N = r->N;
    const int N_active = (r->N_active==SIZE_MAX || r->testparticle_type==1)?r->N:r->N_active;

    const int idxA = 0; // star A assumed to be first particle
    const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary"));

    if (idxB == -1){
        // no binary companion, just shift to DH coordinates
        reb_integrator_trace_dh_to_inertial(r);
        return;
    }

    // ------------------ sums from WB planet coordinates ------------------
    double M = 0.0; // Σ m_i
    struct reb_vec3d sum_mX = {0.0, 0.0, 0.0}; // Σ m_i X_i
    struct reb_vec3d sum_mV = {0.0, 0.0, 0.0}; // Σ m_i V_i

    for (int i = 1; i <= N_active-1; ++i) {
        if (i == idxB) continue; // skip binary companion
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

    // Cache WB star
    const struct reb_vec3d X_B = { particles[idxB].x,  particles[idxB].y,  particles[idxB].z };
    const struct reb_vec3d V_B = { particles[idxB].vx, particles[idxB].vy, particles[idxB].vz };

    // POSITIONS
    struct reb_vec3d sum_mX_over_d = {0.0, 0.0, 0.0};
    sum_mX_over_d.x = sum_mX.x / denomA;
    sum_mX_over_d.y = sum_mX.y / denomA;
    sum_mX_over_d.z = sum_mX.z / denomA;

    // Use com to calculate central object's position.
    // This ignores previous values stored in particles[0].
    // Should not matter unless collisions occurred.
    particles[idxA].x = trace->com_pos.x - (mB * X_B.x + mB * sum_mX_over_d.x + sum_mX.x) / mtot;
    particles[idxA].y = trace->com_pos.y - (mB * X_B.y + mB * sum_mX_over_d.y + sum_mX.y) / mtot;
    particles[idxA].z = trace->com_pos.z - (mB * X_B.z + mB * sum_mX_over_d.z + sum_mX.z) / mtot;

    // Planets: x_i = X_i + x_A
    // Loop over all N to include test particles
    for (int i=1;i<N;i++) {
        if (i == idxB) continue; // skip binary companion
        particles[i].x += particles[idxA].x; // particles[i].x currently X_i
        particles[i].y += particles[idxA].y;
        particles[i].z += particles[idxA].z;
    }

    // Star B: x_B = X_B + x_A + sum_mX/denomA
    particles[idxB].x = X_B.x + particles[idxA].x + sum_mX_over_d.x; 
    particles[idxB].y = X_B.y + particles[idxA].y + sum_mX_over_d.y; 
    particles[idxB].z = X_B.z + particles[idxA].z + sum_mX_over_d.z;

    // Common additive term for all planets: W = V_A - coefB * V_B
    const struct reb_vec3d W = {
        trace->com_vel.x - (mB / denomA) * V_B.x,
        trace->com_vel.y - (mB / denomA) * V_B.y,
        trace->com_vel.z - (mB / denomA) * V_B.z
    };

    // v_A
    particles[idxA].vx = W.x - sum_mV.x / mA;
    particles[idxA].vy = W.y - sum_mV.y / mA;
    particles[idxA].vz = W.z - sum_mV.z / mA;

    // v_B
    particles[idxB].vx = V_B.x + trace->com_vel.x;
    particles[idxB].vy = V_B.y + trace->com_vel.y;
    particles[idxB].vz = V_B.z + trace->com_vel.z;

    // Planets: v_i = V_i + add_planet
    for (int i = 1; i < N; ++i) {
        if (i == idxB) continue;
        particles[i].vx += W.x; // particles[i].vx currently V_i
        particles[i].vy += W.y;
        particles[i].vz += W.z;
    }
}

static void reb_integrator_trace_calculate_acceleration_mode_interaction_wb(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const size_t N_active = ((r->N_active==SIZE_MAX)?N:r->N_active);
    const int _testparticle_type   = r->testparticle_type;

    const int idxA = 0;
    const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary"));
    const int has_binary = (idxB != -1);
    const double mA = particles[idxA].m;
    const double mB = has_binary ? particles[idxB].m : 0.0;

    // little s, equation (13)
    struct reb_vec3d s = {0};
    double mtotA = mA;
    for (int i=0; i<N; i++){
        particles[i].ax = 0;
        particles[i].ay = 0;
        particles[i].az = 0;

        if (i != idxA && i < N_active){ // only loop over active planets
            if (i==idxB) continue; // skip binary companion
            s.x += particles[i].m * particles[i].x;
            s.y += particles[i].m * particles[i].y;
            s.z += particles[i].m * particles[i].z;
            mtotA += particles[i].m;
        }
    }

    s.x /= mtotA;
    s.y /= mtotA;
    s.z /= mtotA;

    // Binary-related quantities (only if binary exists)
    double dbx = 0.0, dby = 0.0, dbz = 0.0;
    double _rb = 1.0, _Rb = 1.0;
    double prefact1 = 0.0;

    if (has_binary){
        dbx = particles[idxB].x + s.x;
        dby = particles[idxB].y + s.y;
        dbz = particles[idxB].z + s.z;
        _rb = sqrt(dbx*dbx + dby*dby + dbz*dbz);
        _Rb = sqrt(particles[idxB].x*particles[idxB].x + particles[idxB].y*particles[idxB].y + particles[idxB].z*particles[idxB].z);
        prefact1 = -G*mA*mB/mtotA / (_rb*_rb*_rb);

        // Acceleration of the binary due to star A
        const double prefactorb1 = G*mA / (_Rb*_Rb*_Rb);
        const double prefactorb2 = G*mA / (_rb*_rb*_rb);
        particles[idxB].ax += prefactorb1*particles[idxB].x - prefactorb2*dbx;
        particles[idxB].ay += prefactorb1*particles[idxB].y - prefactorb2*dby;
        particles[idxB].az += prefactorb1*particles[idxB].z - prefactorb2*dbz;
    }

    // Second term in Equation 14 compute C = sum_i m_i * (Xb - Xi + s)/|...|^3
    struct reb_vec3d C = {0};

    // Active-active planet interactions
    for (int i=2; i<N_active;i++){
        if (reb_sigint > 1) return;
        if (i == idxB) continue;
        for (int j=1;j<i;j++){
            if (j == idxB) continue;
            // Pairwise planet interactions
            if (trace->current_Ks[j*r->N+i]) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz);

            // Last term of Equation (14)
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

    // Repurposing this loop. Outer loop is used to add in the binary terms too.
    const int startitestp = MAX(N_active,2);
    for (int i=1; i<N; i++){
        if (reb_sigint > 1) return;
        if (i == idxB) continue; // to avoid double counting the binary

        // All the binary interactions
        // Xb - Xi + Sx
        if (has_binary){
            // Xb - Xi + Sx
            const double dbix = dbx - particles[i].x;
            const double dbiy = dby - particles[i].y;
            const double dbiz = dbz - particles[i].z;
            const double _rbi = sqrt(dbix*dbix + dbiy*dbiy + dbiz*dbiz);

            const double rx = (particles[idxB].x + s.x) - particles[i].x;
            const double ry = (particles[idxB].y + s.y) - particles[i].y;
            const double rz = (particles[idxB].z + s.z) - particles[i].z;
            const double _rs  = sqrt(rx*rx + ry*ry + rz*rz);
            const double invr3 = 1.0/(_rs*_rs*_rs);

            C.x += particles[i].m * rx * invr3;
            C.y += particles[i].m * ry * invr3;
            C.z += particles[i].m * rz * invr3;

            // accelerations on particle i due to binary
            const double prefact3 = G*mB / (_rbi*_rbi*_rbi);
            particles[i].ax += prefact1*dbx + prefact3*dbix;
            particles[i].ay += prefact1*dby + prefact3*dbiy;
            particles[i].az += prefact1*dbz + prefact3*dbiz;

            // acceleration on the binary due to i
            const double prefactorb1i = G*particles[i].m/(_Rb*_Rb*_Rb);
            const double prefactorb2i = G*particles[i].m/(_rbi*_rbi*_rbi);
            particles[idxB].ax += prefactorb1i*particles[idxB].x - prefactorb2i*dbix;
            particles[idxB].ay += prefactorb1i*particles[idxB].y - prefactorb2i*dbiy;
            particles[idxB].az += prefactorb1i*particles[idxB].z - prefactorb2i*dbiz;
        }

        // Inner loop for active-test particle interactions
        if (i >= startitestp){
            for (int j=1; j<N_active; j++){
                if (j == idxB) continue;
                if (trace->current_Ks[j*r->N+i]) continue;
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
    }

    if (has_binary){
        const double common = -G * mB / mtotA;
        for (int k=1; k<N; k++){
            if (k == idxB) continue;
            particles[k].ax += common * C.x;
            particles[k].ay += common * C.y;
            particles[k].az += common * C.z;
        }
    }

}

static void reb_integrator_trace_calculate_acceleration_mode_interaction(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
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
            if (trace->current_Ks[j*N+i]) continue;
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
            if (trace->current_Ks[j*N+i]) continue;
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
            if (trace->current_Ks[j*N+i]) continue;
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
                if (trace->current_Ks[j*N+i]) continue;
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
        if(r->N>trace->N_allocated_additional_forces){
            trace->particles_backup_additional_forces = realloc(trace->particles_backup_additional_forces, r->N*sizeof(struct reb_particle));
            trace->N_allocated_additional_forces = r->N;
        }
        memcpy(trace->particles_backup_additional_forces,r->particles,r->N*sizeof(struct reb_particle));
        if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_wb_to_inertial(r);
        else reb_integrator_trace_dh_to_inertial(r);
        r->additional_forces(r);
        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle* restrict const backup = trace->particles_backup_additional_forces;
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

//static void reb_integrator_trace_calculate_acceleration_mode_kepler_wb(struct reb_simulation* r){
    
//}

static void reb_integrator_trace_calculate_acceleration_mode_kepler(struct reb_simulation* r){
    // Kepler Step
    struct reb_particle* const particles = r->particles;
    const size_t N = r->N;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    const double m0 = r->particles[0].m;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const size_t encounter_N = trace->encounter_N;
    const size_t encounter_N_active = trace->encounter_N_active;
    size_t* map = trace->encounter_map;
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
                if (!trace->current_Ks[mj*N+mi]) continue;
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
            if (!trace->current_Ks[mj*N+mi]) continue;
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
            if (!trace->current_Ks[mj*N+mi]) continue;
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
                if (!trace->current_Ks[mj*N+mi]) continue;
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
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const size_t N = r->N;
    trace->mode = REB_INTEGRATOR_TRACE_MODE_INTERACTION;
    if (trace->coordinates==REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY){
        reb_integrator_trace_calculate_acceleration_mode_interaction_wb(r);
    }else{
        reb_integrator_trace_calculate_acceleration_mode_interaction(r);
    }

    for (size_t i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_trace_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const int current_C = trace->current_C;
    if (current_C) return; // No jump step for pericenter approaches

    const size_t N_active = r->N_active==SIZE_MAX?r->N:r->N_active;

    // If TP type 1, use r->N. Else, use N_active.
    const size_t N = r->testparticle_type==0 ? N_active: r->N;

    const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary"));

    double px=0., py=0., pz=0.;
    for (size_t i=1;i<N;i++){
        if (trace->coordinates==REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY && i==idxB) continue; // skip star B
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m;
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px *= dt/r->particles[0].m;
    py *= dt/r->particles[0].m;
    pz *= dt/r->particles[0].m;

    const size_t N_all = r->N;
    for (size_t i=1;i<N_all;i++){
        if (trace->coordinates==REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY && i==idxB) continue; // skip star B
        particles[i].x += px;
        particles[i].y += py;
        particles[i].z += pz;
    }
}

void reb_integrator_trace_com_step(struct reb_simulation* const r, double dt){
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    trace->com_pos.x += dt*trace->com_vel.x;
    trace->com_pos.y += dt*trace->com_vel.y;
    trace->com_pos.z += dt*trace->com_vel.z;
}

void reb_integrator_trace_whfast_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const size_t N = r->N;
    const int N_active = r->N_active==SIZE_MAX?r->N:r->N_active;
    
    const int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary"));
    const int has_binary = (idxB != -1);
    struct reb_particle* p = r->particles;

    double mA = p[0].m;
    double mB = has_binary ? p[idxB].m : 0.0;
    double Mpl = 0.0;
    for (size_t i=1; i<N; ++i){
        if (trace->coordinates==REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY && i==idxB){
            // This is ONLY to integrate the wide binary in WB coordinates
            // Kepler parameter for relative coordinate X_B is G * m_tot
            const double GM = r->G * (mA + Mpl + mB);

            // I believe we need to convert from V_B equiv P_B / m_B to real velocities here, just for the Kepler solver.
            // Works! note to self get this sanity checked...

            const double scale = ((mA + mB + Mpl) / (mA + Mpl));
            p[i].vx *= scale;  p[i].vy *= scale;  p[i].vz *= scale;
            reb_integrator_whfast_kepler_solver(&p[i],GM,dt,NULL);
            p[i].vx /= scale;  p[i].vy /= scale;  p[i].vz /= scale;
        }else{
            // normal Keplerian orbit
            reb_integrator_whfast_kepler_solver(&p[i],r->G*p[0].m,dt,NULL);
            Mpl += p[i].m;
        }
    }
}

void reb_integrator_trace_update_particles(struct reb_simulation* r, const double* y){
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    size_t N = trace->encounter_N;
    size_t* map = trace->encounter_map;

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
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    // TRACE always needs this to ensure the right Hamiltonian is evolved
    reb_integrator_trace_update_particles(r, y);
    //if (trace->coordinates==REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY){
    //    reb_integrator_trace_calculate_acceleration_mode_kepler_wb(r);
    //}else{
    reb_integrator_trace_calculate_acceleration_mode_kepler(r);
    //}

    double px=0., py=0., pz=0.;
    size_t* map = trace->encounter_map;
    size_t N = trace->encounter_N;

    if (map==NULL){
        reb_simulation_error(r, "Cannot access TRACE map from BS.");
        return;
    }

    // Kepler Step
    // This is only for pericenter approach
    if (trace->current_C){
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
    struct reb_integrator_trace_state* const trace = r->integrator.state;

    if (trace->encounter_N < 2){
        // No close encounters, skip
        return;
    }

    size_t i_enc = 0;
    const size_t N_active = r->N_active==SIZE_MAX ? r->N : r->N_active;
    trace->encounter_N_active = 0;
    for (size_t i=0; i<r->N; i++){
        if(trace->encounter_map[i]){
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = trace->particles_backup_kepler[i]; // Coordinates before WHFast step, overwrite particles with close encounters
            trace->encounter_map[i_enc] = i;
            i_enc++;
            if (i<N_active){
                trace->encounter_N_active++;
                if (trace->tponly_encounter){
                    trace->particles_backup_kepler[i] = tmp;         // Make copy of particles after the kepler step.
                                                                     // used to restore the massive objects' states in the case
                                                                     // of only massless test-particle encounters
                }
            }
        }
    }

    trace->mode = REB_INTEGRATOR_TRACE_MODE_KEPLER;
    r->map = trace->encounter_map; // for collision search
    r->N_map = trace->encounter_N;
    r->gravity = REB_GRAVITY_CUSTOM;
    //if (trace->coordinates==REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY){
        //r->gravity_custom = reb_integrator_trace_calculate_acceleration_mode_kepler_wb;
    //}else{
    r->gravity_custom = reb_integrator_trace_calculate_acceleration_mode_kepler;
    //}

    // Only Partial BS uses this step 
    if (trace->peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS || !trace->current_C){
        // run
        const double old_dt = r->dt;
        const double old_t = r->t;
        const double t_needed = r->t + dt;
        struct reb_integrator_bs_state* bs = reb_integrator_bs.create();

        // Temporarily remove all odes for BS step
        struct reb_ode** odes_backup = r->odes;
        size_t N_allocated_odes_backup = r->N_allocated_odes;
        size_t N_odes_backup = r->N_odes;
        r->odes = NULL;
        r->N_allocated_odes = 0;
        r->N_odes = 0;

        // Temporarily add new nbody ode for BS step
        struct reb_ode* nbody_ode = NULL;

        // TODO: Support backwards integrations
        while(r->t < t_needed && fabs(dt/old_dt)>1e-14 && r->status<=0){
            if (!nbody_ode || nbody_ode->length != trace->encounter_N*3*2){
                // (re)create the ODE
                reb_ode_free(nbody_ode);
                nbody_ode = reb_ode_create(r, trace->encounter_N*3*2);
                nbody_ode->derivatives = reb_integrator_trace_nbody_derivatives;
                nbody_ode->needs_nbody = 0;
                bs->first_or_last_step = 1;
            }

            double* y = nbody_ode->y;

            // In case of overshoot
            if (r->t + dt >  t_needed){
                dt = t_needed - r->t;
            }

            struct reb_particle star = r->particles[0]; // backup velocity
            r->particles[0].vx = 0; // star does not move in dh
            r->particles[0].vy = 0;
            r->particles[0].vz = 0;

            for (size_t i=0; i<trace->encounter_N; i++){
                const size_t mi = trace->encounter_map[i];
                const struct reb_particle p = r->particles[mi];
                y[i*6+0] = p.x;
                y[i*6+1] = p.y;
                y[i*6+2] = p.z;
                y[i*6+3] = p.vx;
                y[i*6+4] = p.vy;
                y[i*6+5] = p.vz;
            }

            int success = reb_integrator_bs_step_odes(r, bs, dt);
            if (success){
                r->t += dt;
            }
            dt = bs->dt_proposed;
            reb_integrator_trace_update_particles(r, nbody_ode->y);

            r->particles[0].vx = star.vx; // restore every timestep for collisions
            r->particles[0].vy = star.vy;
            r->particles[0].vz = star.vz;

            if (success){
                // Only do a collision search for accepted steps.
                reb_collision_search(r);
                if (r->N_collisions) trace->force_accept = 1;
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
        if(trace->tponly_encounter){
            for (size_t i=1; i < trace->encounter_N_active; i++){
                size_t mi = trace->encounter_map[i];
                r->particles[mi] = trace->particles_backup_kepler[mi];
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
        reb_integrator_bs.free(bs);
    }
    r->map = NULL; // for collision search
    r->N_map = 0;
}

void reb_integrator_trace_kepler_step(struct reb_simulation* const r, const double _dt){
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    memcpy(trace->particles_backup_kepler,r->particles,r->N*sizeof(struct reb_particle));
    reb_integrator_trace_whfast_step(r, _dt);
    reb_integrator_trace_bs_step(r, _dt);
}


void reb_integrator_trace_pre_ts_check(struct reb_simulation* const r){
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const size_t N = r->N;
    const size_t Nactive = r->N_active==SIZE_MAX?r->N:r->N_active;
    int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary"));
    int (*_switch) (struct reb_simulation* const r, const size_t i, const size_t j) = trace->S ? trace->S : reb_integrator_trace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const size_t j) = trace->S_peri ? trace->S_peri : reb_integrator_trace_switch_peri_default;

    // Clear encounter map
    for (size_t i=1; i<r->N; i++){
        trace->encounter_map[i] = 0;
    }
    trace->encounter_map[0] = 1;
    trace->encounter_N = 1;

    // Reset encounter triggers.
    trace->current_C = 0;

    for (size_t i = 0; i < N; i++){
        for (size_t j = i + 1; j < N; j++){
            trace->current_Ks[i*N+j] = 0;
        }
    }

    if (r->testparticle_type == 1){
        trace->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        trace->tponly_encounter = 1;
    }

    // Check for pericenter CE
    // in WB we do the pericenter checks in inertial
    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_wb_to_inertial(r);
    for (size_t j = 1; j < Nactive; j++){
        if (_switch_peri(r, j)){
            trace->current_C = 1;
            if (trace->peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS || trace->peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15){
                // Everything will be integrated with BS/IAS15. No need to check any further.
                if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_inertial_to_wb(r);
                return;
            }
            if (j < Nactive){ // Two massive particles have a close encounter
                trace->tponly_encounter = 0;
                break; // No need to check other particles
            }
        }
    }
    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_inertial_to_wb(r);

    if (trace->current_C){
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        trace->encounter_N = N;
        for (size_t i = 1; i < N; i++){
            trace->encounter_map[i] = 1; //  trigger encounter
        }

    }

    // Body-body
    // there cannot be TP-TP CEs
    for (size_t i = 0; i < Nactive; i++){ // Check central body, for collisions
        for (size_t j = i + 1; j < N; j++){
            if ((i == idxB || j == idxB) && trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) continue; // in WB coordinates star B does not have close encounters, for now
            if (_switch(r, i, j)){
                trace->current_Ks[i*N+j] = 1;
                if (trace->encounter_map[i] == 0){
                    trace->encounter_map[i] = 1; // trigger encounter
                    trace->encounter_N++;
                }
                if (trace->encounter_map[j] == 0){
                    trace->encounter_map[j] = 1; // trigger encounter
                    trace->encounter_N++;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    trace->tponly_encounter = 0;
                }
            }
        }
    }
    memcpy(trace->encounter_map_backup, trace->encounter_map, N*sizeof(size_t));
}

double reb_integrator_trace_post_ts_check(struct reb_simulation* const r){
    // This function returns 1 if any new encounters occurred.
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    const size_t N = r->N;
    const size_t Nactive = r->N_active==SIZE_MAX?r->N:r->N_active;
    int idxB = reb_simulation_particle_index(reb_simulation_get_particle_by_name(r, "widebinary")); // only used for WB
    int (*_switch) (struct reb_simulation* const r, const size_t i, const size_t j) = trace->S ? trace->S : reb_integrator_trace_switch_default;
    int (*_switch_peri) (struct reb_simulation* const r, const size_t j) = trace->S_peri ? trace->S_peri : reb_integrator_trace_switch_peri_default;
    size_t new_close_encounter = 0; // New CEs

    // Set this from pre-ts encounter map. I don't think we need to reset encounter_N here.
    memcpy(trace->encounter_map, trace->encounter_map_backup, N*sizeof(size_t));

    // in WB we do the pericenter checks in inertial
    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_wb_to_inertial(r);
    if (!trace->current_C){
        // Check for pericenter CE if not already triggered from pre-timestep.
        for (size_t j = 1; j < Nactive; j++){
            
            // in WB coordinates this is checked against the WB itself
            if (j == idxB && trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) continue;
            
            if (_switch_peri(r, j)){
                trace->current_C = 1;
                new_close_encounter = 1;
                if (trace->peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS || trace->peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15){
                    // Everything will be integrated with BS/IAS15. No need to check any further.
                    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_inertial_to_wb(r);
                    return new_close_encounter;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    trace->tponly_encounter = 0;
                    break; // No need to check other particles
                }
            }
        }
    }
    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_inertial_to_wb(r);

    if (trace->current_C){
        // Pericenter close encounter detected. We integrate the entire simulation with BS
        trace->encounter_N = N;
        for (size_t i = 0; i < N; i++){
            trace->encounter_map[i] = 1; // trigger encounter
        }
    }


    // Body-body
    // there cannot be TP-TP CEs
    for (size_t i = 0; i < Nactive; i++){ // Do not check for central body anymore
        for (size_t j = i + 1; j < N; j++){
            if ((i == idxB || j == idxB) && trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) continue; // in WB coordinates star B does not have close encounters, for now
            if (_switch(r, i, j)){
                if (trace->current_Ks[i*N+j] == 0){
                    new_close_encounter = 1;
                }
                trace->current_Ks[i*N+j] = 1;
                if (trace->encounter_map[i] == 0){
                    trace->encounter_map[i] = 1; // trigger encounter
                    trace->encounter_N++;
                }
                if (trace->encounter_map[j] == 0){
                    trace->encounter_map[j] = 1; // trigger encounter
                    trace->encounter_N++;
                }

                if (j < Nactive){ // Two massive particles have a close encounter
                    trace->tponly_encounter = 0;
                }
            }
        }
    }

    return new_close_encounter;
}

static void reb_integrator_trace_step_try(struct reb_simulation* const r){
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    if (trace->current_C == 0 || trace->peri_mode == REB_INTEGRATOR_TRACE_PERIMODE_PARTIAL_BS){
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
        trace->mode = REB_INTEGRATOR_TRACE_MODE_FULL;
        if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_wb_to_inertial(r);
        else reb_integrator_trace_dh_to_inertial(r);
        switch (trace->peri_mode){
            case REB_INTEGRATOR_TRACE_PERIMODE_FULL_IAS15:
                {
                    struct reb_integrator_ias15_state* ias15 = reb_integrator_ias15.create();
                    while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 && r->status<=0){
                        reb_integrator_ias15.step(r, ias15);
                        if (r->t+r->dt >  t_needed){
                            r->dt = t_needed-r->t;
                        }
                        reb_collision_search(r);
                        if (r->N_collisions) trace->force_accept = 1;
                    }
                    reb_integrator_ias15.free(ias15);
                }
                break;
            case REB_INTEGRATOR_TRACE_PERIMODE_FULL_BS:
                {
                    struct reb_integrator_bs_state* bs = reb_integrator_bs.create();
                    struct reb_ode* nbody_ode = NULL;

                    double* y;
                    while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 && r->status<=0){
                        if (!nbody_ode || nbody_ode->length != 6*r->N){
                            // (re)create the ODE
                            reb_ode_free(nbody_ode);
                            nbody_ode = reb_ode_create(r, 6*r->N);
                            nbody_ode->derivatives = reb_integrator_bs_nbody_derivatives;
                            nbody_ode->needs_nbody = 0;
                            bs->first_or_last_step = 1;
                            y = nbody_ode->y;
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

                        int success = reb_integrator_bs_step_odes(r, bs, r->dt);
                        if (success){
                            r->t += r->dt;
                        }
                        r->dt = bs->dt_proposed;
                        if (r->t+r->dt >  t_needed){
                            r->dt = t_needed-r->t;
                        }

                        reb_integrator_bs_update_particles(r, nbody_ode->y);

                        if (success){
                            // Only do a collision search for accepted steps.
                            reb_collision_search(r);
                            if (r->N_collisions) trace->force_accept = 1;
                        }
                    }
                    reb_ode_free(nbody_ode);
                    reb_integrator_bs.free(bs);
                }
                break;
            default:
                reb_simulation_error(r,"Unsupported peri_mode encountered\n");
                break;
        }
        r->t = old_t; // final time will be set later
        r->dt = old_dt;
        
        if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_inertial_to_wb(r);
        else reb_integrator_trace_inertial_to_dh(r);
    }
}

void reb_integrator_trace_did_add_particle(struct reb_simulation* r){
    // TRACE can add particles mid-timestep now
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    if (trace->mode==REB_INTEGRATOR_TRACE_MODE_KEPLER){
        const size_t old_N = r->N-1;
        if (trace->N_allocated < r->N){
            trace->current_Ks    = realloc(trace->current_Ks, sizeof(int)*r->N*r->N);
            trace->particles_backup = realloc(trace->particles_backup, sizeof(struct reb_particle)*r->N);
            trace->particles_backup_kepler = realloc(trace->particles_backup_kepler, sizeof(struct reb_particle)*r->N);
            trace->current_Ks    = realloc(trace->current_Ks, sizeof(int)*r->N*r->N);
            trace->encounter_map = realloc(trace->encounter_map, sizeof(size_t)*r->N);
            r->map = trace->encounter_map;
            trace->encounter_map_backup = realloc(trace->encounter_map_backup, sizeof(size_t)*r->N);
            trace->N_allocated   = r->N;
        }

        // First reshuffle existing Ks
        size_t i = old_N;
        while (i --> 0){
            size_t j = old_N;
            while (j --> 0){
                trace->current_Ks[i*old_N+j+i] = trace->current_Ks[i*old_N+j];
            }
        }

        // add in new particle, we want it to interact with all currently interacting particles
        // exclude star
        for (size_t i = 1; i < trace->encounter_N; i++){
            trace->current_Ks[trace->encounter_map[i]*r->N+old_N] = 1;
        }

        trace->encounter_map[trace->encounter_N] = old_N;
        trace->encounter_N++;
        r->N_map++;

        if (r->N_active==SIZE_MAX){ 
            // If global N_active is not set, then all particles are active, so the new one as well.
            // Otherwise, assume we're adding non active particle. 
            trace->encounter_N_active++;
        }

    }
}

void reb_integrator_trace_will_remove_particle(struct reb_simulation* r, size_t index){
    struct reb_integrator_trace_state* const trace = r->integrator.state;
    // TODO: Logic unclear here. A BS reset might still need to be reimplemented here.
    //reb_integrator_bs_reset(r);
    if (trace->mode==REB_INTEGRATOR_TRACE_MODE_KEPLER){
        // Only removed mid-timestep if collision - BS Step!
        int after_to_be_removed_particle = 0;
        size_t encounter_index = SIZE_MAX;
        for (size_t i=0;i<trace->encounter_N;i++){
            if (after_to_be_removed_particle == 1){
                trace->encounter_map[i-1] = trace->encounter_map[i] - 1;
            }
            if (trace->encounter_map[i]==index){
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
                trace->current_Ks[i*new_N+j] = trace->current_Ks[i*new_N+j+counter];
            }
        }
        if (encounter_index<trace->encounter_N_active){
            trace->encounter_N_active--;
        }
        trace->encounter_N--;
        r->N_map--;
    }
}


void reb_integrator_trace_step(struct reb_simulation* r, void* state){
    // Do memory management and consistency checks
    struct reb_integrator_trace_state* const trace = state;
    const size_t N = r->N;

    if (r->N_var){
        reb_simulation_warning(r,"TRACE does not work with variational equations.");
    }

    if (trace->N_allocated<N){
        // These arrays are only used within one timestep.
        // Can be recreated without loosing bit-wise reproducibility.
        trace->particles_backup       = realloc(trace->particles_backup,sizeof(struct reb_particle)*N);
        trace->particles_backup_kepler   = realloc(trace->particles_backup_kepler,sizeof(struct reb_particle)*N);
        trace->current_Ks             = realloc(trace->current_Ks,sizeof(int)*N*N);
        trace->encounter_map          = realloc(trace->encounter_map,sizeof(size_t)*N);
        trace->encounter_map_backup   = realloc(trace->encounter_map_backup,sizeof(size_t)*N);
        trace->N_allocated = N;
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

    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_inertial_to_wb(r);
    else reb_integrator_trace_inertial_to_dh(r);

    // Create copy of all particle to allow for the step to be rejected.
    memcpy(trace->particles_backup, r->particles, N*sizeof(struct reb_particle));

    // This will be set to 1 if a collision occurred.
    trace->force_accept = 0;

    // Check if there are any close encounters
    reb_integrator_trace_pre_ts_check(r);

    // Attempt one step. 
    reb_integrator_trace_step_try(r);

    // We always accept the step if a collision occurred as it is impossible to undo the collision.
    if (!trace->force_accept){
        // We check again for close encounters to ensure time reversibility. 
        if (reb_integrator_trace_post_ts_check(r)){
            // New encounters were found. Will reject the step.
            // Revert particles to the beginning of the step.
            memcpy(r->particles, trace->particles_backup, N*sizeof(struct reb_particle));

            // Do step again
            reb_integrator_trace_step_try(r);
        }
    }
    if (trace->coordinates == REB_INTEGRATOR_TRACE_COORDINATES_WIDEBINARY) reb_integrator_trace_wb_to_inertial(r);
    else reb_integrator_trace_dh_to_inertial(r);

    r->t+=r->dt;
    r->dt_last_done = r->dt;
    r->N_targets = 1; // Only search for collisions with star after complete timestep.
}

