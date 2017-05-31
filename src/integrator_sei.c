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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "rebound.h"
#include "particle.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator.h"
#include "integrator_sei.h"


static void operator_H012(double dt, const struct reb_simulation_integrator_sei ri_sei, struct reb_particle* p);
static void operator_phi1(double dt, struct reb_particle* p);


void reb_integrator_sei_init(struct reb_simulation* const r){
    /**
     * Pre-calculates sin() and tan() needed for SEI. 
     */
	if (r->ri_sei.OMEGAZ==-1){
		r->ri_sei.OMEGAZ=r->ri_sei.OMEGA;
	}
    r->ri_sei.sindt = sin(r->ri_sei.OMEGA*(-r->dt/2.));
    r->ri_sei.tandt = tan(r->ri_sei.OMEGA*(-r->dt/4.));
    r->ri_sei.sindtz = sin(r->ri_sei.OMEGAZ*(-r->dt/2.));
    r->ri_sei.tandtz = tan(r->ri_sei.OMEGAZ*(-r->dt/4.));
    r->ri_sei.lastdt = r->dt;
}

void reb_integrator_sei_part1(struct reb_simulation* const r){
    r->gravity_ignore_terms = 0;
	const int N = r->N;
	struct reb_particle* const particles = r->particles;
	if (r->ri_sei.OMEGAZ==-1){
		r->ri_sei.OMEGAZ=r->ri_sei.OMEGA;
	}
	if (r->ri_sei.lastdt!=r->dt){
        reb_integrator_sei_init(r);
	}
	const struct reb_simulation_integrator_sei ri_sei = r->ri_sei;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		operator_H012(r->dt, ri_sei, &(particles[i]));
	}
	r->t+=r->dt/2.;
}

void reb_integrator_sei_part2(struct reb_simulation* r){
	const int N = r->N;
	struct reb_particle* const particles = r->particles;
	const struct reb_simulation_integrator_sei ri_sei = r->ri_sei;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		operator_phi1(r->dt, &(particles[i]));
		operator_H012(r->dt, ri_sei, &(particles[i]));
	}
	r->t+=r->dt/2.;
	r->dt_last_done = r->dt;
}

void reb_integrator_sei_synchronize(struct reb_simulation* r){
	// Do nothing.
}

void reb_integrator_sei_reset(struct reb_simulation* r){
	r->ri_sei.lastdt = 0;	
}

/**
 * @brief This function evolves a particle under the unperturbed
 * Hamiltonian H0 exactly up to machine precission.
 * @param p reb_particle to evolve.
 * @param dt Timestep
 * @param ri_sei Integrator struct
 */
static void operator_H012(double dt, const struct reb_simulation_integrator_sei ri_sei, struct reb_particle* p){
		
	// Integrate vertical motion
	const double zx = p->z * ri_sei.OMEGAZ;
	const double zy = p->vz;
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double zt1 =  zx - ri_sei.tandtz*zy;			
	const double zyt =  ri_sei.sindtz*zt1 + zy;
	const double zxt =  zt1 - ri_sei.tandtz*zyt;	
	p->z  = zxt/ri_sei.OMEGAZ;
	p->vz = zyt;

	// Integrate motion in xy directions
	const double aO = 2.*p->vy + 4.*p->x*ri_sei.OMEGA;	// Center of epicyclic motion
	const double bO = p->y*ri_sei.OMEGA - 2.*p->vx;	

	const double ys = (p->y*ri_sei.OMEGA-bO)/2.; 		// Epicycle vector
	const double xs = (p->x*ri_sei.OMEGA-aO); 
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double xst1 =  xs - ri_sei.tandt*ys;			
	const double yst  =  ri_sei.sindt*xst1 + ys;
	const double xst  =  xst1 - ri_sei.tandt*yst;	

	p->x  = (xst+aO)    /ri_sei.OMEGA;			
	p->y  = (yst*2.+bO) /ri_sei.OMEGA - 3./4.*aO*dt;	
	p->vx = yst;
	p->vy = -xst*2. -3./2.*aO;
}

/**
 * @brief This function applies the acceleration due to the PHI1 term.
 * @details It is only exact if the forces are velocity independet (i.e. gravity).
 * If the forces are velocity dependent, it breaks the symmetry of the scheme,
 * making it firsr-order and non-symplectic. As long as these forces are small,
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

