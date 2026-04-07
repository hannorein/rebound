/**
 * @file 	integrator_sei.c
 * @brief 	Symplectic Epicycle Integrator (SEI).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the Symplectic Epicycle Integrator 
 * (SEI). The integrator is described in detail in Rein & Tremaine 2011. 
 * It solves epicyclic motion exactly and is therefore exact up to machine
 * precision in the limit of no perturbing forces. When perturbing-forces
 * are of order eps, then the error of the scheme is O(eps dt^3). It also
 * makes use of two shear operators instead of a rotation to minimize 
 * systematic numerical round-off errors.
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
#include "rebound.h"
#include "rebound_internal.h"
#include <math.h>
#include <stddef.h>
#include <assert.h>
#include "particle.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_sei.h"
#include "binarydata.h"

// Symplectic Epicycle Integrator SEI (Rein & Tremaine 2011)
struct reb_integrator_sei_state {
    double lastdt;      // Cached sin(), tan() for this value of dt.
    double sindt;       // Cached sin() 
    double tandt;       // Cached tan() 
    double sindtz;      // Cached sin(), z axis
    double tandtz;      // Cached tan(), z axis
};

const struct reb_binarydata_field_descriptor reb_integrator_sei_field_descriptor_list[] = {
    { 56, REB_DOUBLE,       "sei: lastdt",                offsetof(struct reb_integrator_sei_state, lastdt), 0, 0},
    { 57, REB_DOUBLE,       "sei: sindt",                 offsetof(struct reb_integrator_sei_state, sindt), 0, 0},
    { 58, REB_DOUBLE,       "sei: tandt",                 offsetof(struct reb_integrator_sei_state, tandt), 0, 0},
    { 59, REB_DOUBLE,       "sei: sindtz",                offsetof(struct reb_integrator_sei_state, sindtz), 0, 0},
    { 60, REB_DOUBLE,       "sei: tandtz",                offsetof(struct reb_integrator_sei_state, tandtz), 0, 0},
    { 0 }, // Null terminated list
};

const struct reb_integrator reb_integrator_sei = {
    .id = 2,
    .step = reb_integrator_sei_step,
    .reset = reb_integrator_sei_reset,
    .field_descriptor_list = reb_integrator_sei_field_descriptor_list,
};

static void operator_H012(double dt, const struct reb_integrator_sei_state sei, struct reb_particle* p, double OMEGA, double OMEGAZ);
static void operator_phi1(double dt, struct reb_particle* p);


void reb_integrator_sei_alloc(struct reb_simulation* const r){
    assert(!r->integrator_data);
    r->integrator_data = calloc(sizeof(struct reb_integrator_sei_state),1);
}

void reb_integrator_sei_step(struct reb_simulation* const r){
    struct reb_integrator_sei_state* sei = r->integrator_data;
    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    if (sei->lastdt!=r->dt){
        // Pre-calculates sin() and tan() needed for SEI. 
        if (r->OMEGAZ==-1){
            r->OMEGAZ=r->OMEGA;
        }
        sei->sindt = sin(r->OMEGA*(-r->dt/2.));
        sei->tandt = tan(r->OMEGA*(-r->dt/4.));
        sei->sindtz = sin(r->OMEGAZ*(-r->dt/2.));
        sei->tandtz = tan(r->OMEGAZ*(-r->dt/4.));
        sei->lastdt = r->dt;
    }
#pragma omp parallel for schedule(guided)
    for (int i=0;i<N;i++){
        operator_H012(r->dt, *sei, &(particles[i]),r->OMEGA,r->OMEGAZ);
    }
    r->t+=r->dt/2.;

    reb_simulation_update_acceleration(r);

#pragma omp parallel for schedule(guided)
    for (int i=0;i<N;i++){
        operator_phi1(r->dt, &(particles[i]));
        operator_H012(r->dt, *sei, &(particles[i]),r->OMEGA,r->OMEGAZ);
    }
    r->t+=r->dt/2.;
    r->dt_last_done = r->dt;
}

void reb_integrator_sei_reset(struct reb_simulation* r){
    struct reb_integrator_sei_state* sei = r->integrator_data;
    sei->lastdt    = 0;	
}

/**
 * @brief This function evolves a particle under the unperturbed
 * Hamiltonian H0 exactly up to machine precision.
 * @param p reb_particle to evolve.
 * @param dt Timestep
 * @param ri_sei Integrator struct
 */
static void operator_H012(double dt, const struct reb_integrator_sei_state ri_sei, struct reb_particle* p, double OMEGA, double OMEGAZ){

    // Integrate vertical motion
    const double zx = p->z * OMEGAZ;
    const double zy = p->vz;

    // Rotation implemented as 3 shear operators
    // to avoid round-off errors
    const double zt1 =  zx - ri_sei.tandtz*zy;			
    const double zyt =  ri_sei.sindtz*zt1 + zy;
    const double zxt =  zt1 - ri_sei.tandtz*zyt;	
    p->z  = zxt/OMEGAZ;
    p->vz = zyt;

    // Integrate motion in xy directions
    const double aO = 2.*p->vy + 4.*p->x*OMEGA;	// Center of epicyclic motion
    const double bO = p->y*OMEGA - 2.*p->vx;	

    const double ys = (p->y*OMEGA-bO)/2.; 		// Epicycle vector
    const double xs = (p->x*OMEGA-aO); 

    // Rotation implemented as 3 shear operators
    // to avoid round-off errors
    const double xst1 =  xs - ri_sei.tandt*ys;			
    const double yst  =  ri_sei.sindt*xst1 + ys;
    const double xst  =  xst1 - ri_sei.tandt*yst;	

    p->x  = (xst+aO)    /OMEGA;			
    p->y  = (yst*2.+bO) /OMEGA - 3./4.*aO*dt;	
    p->vx = yst;
    p->vy = -xst*2. -3./2.*aO;
}

/**
 * @brief This function applies the acceleration due to the PHI1 term.
 * @details It is only exact if the forces are velocity independent (i.e. gravity).
 * If the forces are velocity dependent, it breaks the symmetry of the scheme,
 * making it first-order and non-symplectic. As long as these forces are small,
 * this should not be visible. However, it is worth keeping in mind. 
 * @param p reb_particle to evolve.
 * @param dt Timestep
 */
static void operator_phi1(double dt, struct reb_particle* p){
    // The force used here is for test cases 2 and 3 
    // in Rein & Tremaine 2011. 
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->vz += p->az * dt;
}

