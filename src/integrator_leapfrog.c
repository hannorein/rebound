/**
 * integrator_leapfrog.c: The standard Leap Frog integator and higher order generalizations
 *
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

#include "rebound.h"
#include "integrator_leapfrog.h"
#include "binarydata.h"

// Leapfrog Integrator (TU splitting)
struct reb_integrator_leapfrog_state {
    unsigned int order;
};
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
    "accurate results when a collision search is turned on."
    ,
    .step = reb_integrator_leapfrog_step,
    .create = reb_integrator_leapfrog_create,
    .free = reb_integrator_leapfrog_free,
    .field_descriptor_list = reb_integrator_leapfrog_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_leapfrog_field_descriptor_list[] = {
    { "Order of the integrator. Default is 2. Other allowed values are 6 and 8.",
        REB_UINT,        "order",          offsetof(struct reb_integrator_leapfrog_state, order), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* reb_integrator_leapfrog_create(){
    struct reb_integrator_leapfrog_state* leapfrog = calloc(sizeof(struct reb_integrator_leapfrog_state),1);
    leapfrog->order = 2;
    return leapfrog;
}

void reb_integrator_leapfrog_free(void* p){
    struct reb_integrator_leapfrog_state* leapfrog = p;
    free(leapfrog);
}


const double reb_integrator_leapfrog_lf4_a = 0.675603595979828817023843904485;
const double reb_integrator_leapfrog_lf6_a[5] = {0.1867, 0.5554970237124784, 0.1294669489134754, -0.843265623387734, 0.9432033015235604};
const double reb_integrator_leapfrog_lf8_a[9] = {0.128865979381443, 0.581514087105251, -0.410175371469850, 0.1851469357165877, -0.4095523434208514, 0.1444059410800120, 0.2783355003936797, 0.3149566839162949, -0.6269948254051343979}; 

static void drift(struct reb_simulation* r, double dt){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
#pragma omp parallel for schedule(guided)
    for (size_t i=0;i<N;i++){
        particles[i].x  += dt * particles[i].vx;
        particles[i].y  += dt * particles[i].vy;
        particles[i].z  += dt * particles[i].vz;
    }
    r->t += dt; // kick step advanced time so that force evaluations are correct.
}

static void kick(struct reb_simulation* r, double dt){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
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
    const double dt = r->dt;
    struct reb_integrator_leapfrog_state* leapfrog = state;
    switch (leapfrog->order){
        case 2:
            drift(r, dt*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt);
            drift(r, dt*0.5);
            break;
        case 4:
            drift(r, dt*reb_integrator_leapfrog_lf4_a);
            reb_simulation_update_acceleration(r);
            kick(r, dt*2.*reb_integrator_leapfrog_lf4_a);
            drift(r, dt*(0.5-reb_integrator_leapfrog_lf4_a));
            reb_simulation_update_acceleration(r);
            kick(r, dt*(1.-4.*reb_integrator_leapfrog_lf4_a));
            drift(r, dt*(0.5-reb_integrator_leapfrog_lf4_a));
            reb_simulation_update_acceleration(r);
            kick(r, dt*2.*reb_integrator_leapfrog_lf4_a);
            drift(r, dt*reb_integrator_leapfrog_lf4_a);
            break;
        case 6:
            drift(r, dt*reb_integrator_leapfrog_lf6_a[0]*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[0]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[0]+reb_integrator_leapfrog_lf6_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[1]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[1]+reb_integrator_leapfrog_lf6_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[2]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[2]+reb_integrator_leapfrog_lf6_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[3]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[3]+reb_integrator_leapfrog_lf6_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[4]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[3]+reb_integrator_leapfrog_lf6_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[3]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[2]+reb_integrator_leapfrog_lf6_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[2]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[1]+reb_integrator_leapfrog_lf6_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[1]);
            drift(r, dt*(reb_integrator_leapfrog_lf6_a[0]+reb_integrator_leapfrog_lf6_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf6_a[0]);
            drift(r, dt*reb_integrator_leapfrog_lf6_a[0]*0.5);
            break; 
        case 8: 
            drift(r, dt*reb_integrator_leapfrog_lf8_a[0]*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[0]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[0]+reb_integrator_leapfrog_lf8_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[1]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[1]+reb_integrator_leapfrog_lf8_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[2]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[2]+reb_integrator_leapfrog_lf8_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[3]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[3]+reb_integrator_leapfrog_lf8_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[4]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[4]+reb_integrator_leapfrog_lf8_a[5])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[5]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[5]+reb_integrator_leapfrog_lf8_a[6])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[6]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[6]+reb_integrator_leapfrog_lf8_a[7])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[7]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[7]+reb_integrator_leapfrog_lf8_a[8])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[8]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[7]+reb_integrator_leapfrog_lf8_a[8])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[7]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[6]+reb_integrator_leapfrog_lf8_a[7])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[6]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[5]+reb_integrator_leapfrog_lf8_a[6])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[5]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[4]+reb_integrator_leapfrog_lf8_a[5])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[4]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[3]+reb_integrator_leapfrog_lf8_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[3]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[2]+reb_integrator_leapfrog_lf8_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[2]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[1]+reb_integrator_leapfrog_lf8_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[1]);
            drift(r, dt*(reb_integrator_leapfrog_lf8_a[0]+reb_integrator_leapfrog_lf8_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, dt*reb_integrator_leapfrog_lf8_a[0]);
            drift(r, dt*reb_integrator_leapfrog_lf8_a[0]*0.5);
            break;
        default:
            reb_simulation_error(r, "Leapfrog order not supported.");
            return;
    }
    r->dt_last_done = dt;
}

