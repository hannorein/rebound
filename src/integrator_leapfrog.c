/**
 * integrator_leapfrog.c: The standard Leap Frog integator, higher order generalizations, and regularization
 *
 * Copyright (c) 2026 Hanno Rein, Stuart Williamson, Shangfei Liu
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
#include "integrator_leapfrog.h"
#include "binarydata.h"
#include <math.h>

void reb_integrator_leapfrog_step(struct reb_simulation* r, void* state);
void* reb_integrator_leapfrog_create();
void reb_integrator_leapfrog_free(void* p);
const struct reb_binarydata_field_descriptor reb_integrator_leapfrog_field_descriptor_list[];

const struct reb_integrator reb_integrator_leapfrog = {
    .documentation = 
    "This is the standard leap frog integrator. It is symplectic. "
    "By default it is second order with one force evaluation per "
    "step. Higher orders of 4, 6, and 8 can be selected as well. "
    "These correspond to the 4th order Yoshida integrator and the "
    "8th order by Blanes & Casa (2016), p91. The higher order methods "
    "have more function evaluations and are therefore slower. Note "
    "that some substeps of the higher order methods move particles "
    "backwards. Therefore higher order methods might not give "
    "accurate results when a collision search is turned on. "
    "A regularized version of leap frog is also available by "
    "setting the regularized flag to 1. This corresponds to the "
    "integrator described by Mikkola and Tanikawa (1999)."
    ,
    .step = reb_integrator_leapfrog_step,
    .create = reb_integrator_leapfrog_create,
    .free = reb_integrator_leapfrog_free,
    .field_descriptor_list = reb_integrator_leapfrog_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_leapfrog_field_descriptor_list[] = {
    { "Order of the integrator. Default is 2. Other allowed values are 6 and 8.",
        REB_UINT,        "order",          offsetof(struct reb_integrator_leapfrog_state, order), 0, 0, 0},
    { "Set this flag to 1 to use the Mikkola and Tanikawa (1999) regularization. Default is 0.",
        REB_UINT,        "regularized",    offsetof(struct reb_integrator_leapfrog_state, regularized), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* reb_integrator_leapfrog_create(){
    struct reb_integrator_leapfrog_state* leapfrog = calloc(sizeof(struct reb_integrator_leapfrog_state),1);
    leapfrog->order = 2;
    leapfrog->regularized = 0;
    leapfrog->E0 = nan("");
    return leapfrog;
}

void reb_integrator_leapfrog_free(void* p){
    struct reb_integrator_leapfrog_state* leapfrog = p;
    free(leapfrog);
}


const double reb_integrator_leapfrog_lf4_a = 0.675603595979828817023843904485;
const double reb_integrator_leapfrog_lf6_a[5] = {0.1867, 0.5554970237124784, 0.1294669489134754, -0.843265623387734, 0.9432033015235604};
const double reb_integrator_leapfrog_lf8_a[9] = {0.128865979381443, 0.581514087105251, -0.410175371469850, 0.1851469357165877, -0.4095523434208514, 0.1444059410800120, 0.2783355003936797, 0.3149566839162949, -0.6269948254051343979}; 

static double potential_energy(struct reb_simulation* r){
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX)?N:r->N_active;
    const struct reb_particle* restrict const particles = r->particles;
    double e_pot = 0.;
    size_t N_interact = (r->testparticle_type==0)?N_active:N;
    for (size_t i=0;i<N_active;i++){
        struct reb_particle pi = particles[i];
        for (size_t j=i+1;j<N_interact;j++){
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;
            e_pot -= r->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    return e_pot;
}

static double kinetic_energy(struct reb_simulation* r){
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX)?N:r->N_active;
    const struct reb_particle* restrict const particles = r->particles;
    double e_kin = 0.;
    size_t N_interact = (r->testparticle_type==0)?N_active:N;
    for (size_t i=0;i<N_interact;i++){
        struct reb_particle pi = particles[i];
        e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
    }
    return e_kin;
}

static double drift(struct reb_simulation* r, struct reb_integrator_leapfrog_state* leapfrog, double dt){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    
    if (leapfrog->regularized){
        dt /= (kinetic_energy(r) - leapfrog->E0); // use ficticious time to determine actual time step
    }

#pragma omp parallel for schedule(guided)
    for (size_t i=0;i<N;i++){
        particles[i].x  += dt * particles[i].vx;
        particles[i].y  += dt * particles[i].vy;
        particles[i].z  += dt * particles[i].vz;
    }
    r->t += dt; // kick step advanced time so that force evaluations are correct.
    return dt;
}

static void kick(struct reb_simulation* r, struct reb_integrator_leapfrog_state* leapfrog, double dt){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;

    if (leapfrog->regularized){
        dt /= -potential_energy(r); // use ficticious time to determine actual time step
    }

#pragma omp parallel for schedule(guided)
    for (size_t i=0;i<N;i++){
        particles[i].vx += dt * particles[i].ax;
        particles[i].vy += dt * particles[i].ay;
        particles[i].vz += dt * particles[i].az;
    }
}

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void reb_integrator_leapfrog_step(struct reb_simulation* r, void* state){
    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    const double dtau = r->dt; // Fictitous timestep for regularized leapfrog, real timestep otherwise
    struct reb_integrator_leapfrog_state* leapfrog = state;

    
    if (leapfrog->regularized){
        if (r->exact_finish_time){
            reb_simulation_error(r, "Regularized leap frog has an adaptive timestep and can not be used with exact_finish_time=1.");
        }
        if (isnan(leapfrog->E0)){
            // If we haven't set the initial energy yet, do it.
            leapfrog->E0 = reb_simulation_energy(r);
        }
    }

    double dt = 0; // real timesteps performed
    switch (leapfrog->order){
        case 2: // 1*2 FE
            dt += drift(r, leapfrog, dtau*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau);
            dt += drift(r, leapfrog, dtau*0.5);
            break;
        case 4: // 3*2 FE
            dt += drift(r, leapfrog, dtau*reb_integrator_leapfrog_lf4_a);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*2.*reb_integrator_leapfrog_lf4_a);
            dt += drift(r, leapfrog, dtau*(0.5-reb_integrator_leapfrog_lf4_a));
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*(1.-4.*reb_integrator_leapfrog_lf4_a));
            dt += drift(r, leapfrog, dtau*(0.5-reb_integrator_leapfrog_lf4_a));
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*2.*reb_integrator_leapfrog_lf4_a);
            dt += drift(r, leapfrog, dtau*reb_integrator_leapfrog_lf4_a);
            break;
        case 6: // 9*2 FE
            dt += drift(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[0]+reb_integrator_leapfrog_lf6_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[1]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[1]+reb_integrator_leapfrog_lf6_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[2]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[2]+reb_integrator_leapfrog_lf6_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[3]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[3]+reb_integrator_leapfrog_lf6_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[4]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[3]+reb_integrator_leapfrog_lf6_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[3]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[2]+reb_integrator_leapfrog_lf6_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[2]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[1]+reb_integrator_leapfrog_lf6_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[1]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf6_a[0]+reb_integrator_leapfrog_lf6_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]);
            dt += drift(r, leapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]*0.5);
            break; 
        case 8: // 17*2 FE
            dt += drift(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[0]+reb_integrator_leapfrog_lf8_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[1]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[1]+reb_integrator_leapfrog_lf8_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[2]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[2]+reb_integrator_leapfrog_lf8_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[3]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[3]+reb_integrator_leapfrog_lf8_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[4]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[4]+reb_integrator_leapfrog_lf8_a[5])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[5]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[5]+reb_integrator_leapfrog_lf8_a[6])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[6]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[6]+reb_integrator_leapfrog_lf8_a[7])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[7]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[7]+reb_integrator_leapfrog_lf8_a[8])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[8]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[7]+reb_integrator_leapfrog_lf8_a[8])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[7]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[6]+reb_integrator_leapfrog_lf8_a[7])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[6]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[5]+reb_integrator_leapfrog_lf8_a[6])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[5]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[4]+reb_integrator_leapfrog_lf8_a[5])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[4]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[3]+reb_integrator_leapfrog_lf8_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[3]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[2]+reb_integrator_leapfrog_lf8_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[2]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[1]+reb_integrator_leapfrog_lf8_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[1]);
            dt += drift(r, leapfrog, dtau*(reb_integrator_leapfrog_lf8_a[0]+reb_integrator_leapfrog_lf8_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]);
            dt += drift(r, leapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]*0.5);
            break;
        default:
            reb_simulation_error(r, "Leap frog order not supported.");
            return;
    }
    r->dt_last_done = dt;
}
