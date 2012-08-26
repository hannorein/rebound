/**
 * @file 	integrator.c
 * @brief 	Symplectic Epicycle Integrator (SEI).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the Symplectic Epicycle Integrator 
 * (SEI). The integrator is described in detail in Rein & Tremaine 2011. 
 * It solves epicyclic motion exactly and is therefore exact up to machine
 * precission in the limit of no perturbing forces. When perturbing-forces
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
#include "particle.h"
#include "gravity.h"
#include "main.h"
#include "boundaries.h"


double OMEGA 	= 1.; 	/**< Epicyclic/orbital frequency. */
double OMEGAZ 	= -1.; 	/**< Epicyclic frequency in vertical direction. */

void operator_H012(struct particle* p);
void operator_phi1(struct particle* p);
// Cache sin() tan() values.
double lastdt=0;	/**< Cached sin(), tan() for this value of dt.*/
double sindt,  tandt;
double sindtz, tandtz;
	
/**
 * This function pre-calculates sin() and tan() needed for SEI. 
 */
void integrator_cache_coefficients(){
	if (lastdt!=dt){
		// Only calculate sin() and tan() if timestep changed
		if (OMEGAZ==-1){
			OMEGAZ=OMEGA;
		}
		sindt = sin(OMEGA*(-dt/2.));
		tandt = tan(OMEGA*(-dt/4.));
		sindtz = sin(OMEGAZ*(-dt/2.));
		tandtz = tan(OMEGAZ*(-dt/4.));
		lastdt = dt;
	}
}

void integrator_part1(){
	integrator_cache_coefficients();
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		operator_H012(&(particles[i]));
	}
	t+=dt/2.;
}

void integrator_part2(){
	integrator_cache_coefficients();
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		operator_phi1(&(particles[i]));
		operator_H012(&(particles[i]));
	}
	t+=dt/2.;
}

/**
 * This function evolves a particle under the unperturbed
 * Hamiltonian H0 exactly up to machine precission.
 * @param p Particle to evolve.
 */
void operator_H012(struct particle* p){
		
	// Integrate vertical motion
	const double zx = p->z * OMEGAZ;
	const double zy = p->vz;
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double zt1 =  zx - tandtz*zy;			
	const double zyt =  sindtz*zt1 + zy;
	const double zxt =  zt1 - tandtz*zyt;	
	p->z  = zxt/OMEGAZ;
	p->vz = zyt;

	// Integrate motion in xy directions
	const double aO = 2.*p->vy + 4.*p->x*OMEGA;	// Center of epicyclic motion
	const double bO = p->y*OMEGA - 2.*p->vx;	

	const double ys = (p->y*OMEGA-bO)/2.; 		// Epicycle vector
	const double xs = (p->x*OMEGA-aO); 
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double xst1 =  xs - tandt*ys;			
	const double yst  =  sindt*xst1 + ys;
	const double xst  =  xst1 - tandt*yst;	

	p->x  = (xst+aO)    /OMEGA;			
	p->y  = (yst*2.+bO) /OMEGA - 3./4.*aO*dt;	
	p->vx = yst;
	p->vy = -xst*2. -3./2.*aO;
}

/**
 * This function applies the acceleration.
 * It is only exact if the forces are velocity independet (i.e. gravity).
 * If the forces are velocity dependent, it breaks the symmetry of the scheme,
 * making it firsr-order and non-symplectic. As long as these forces are small,
 * this should not be visible. However, it is worth keeping in mind. 
 * @param p Particle to evolve.
 */
void operator_phi1(struct particle* p){
	// The force used here is for test cases 2 and 3 
	// in Rein & Tremaine 2011. 
	p->vx += p->ax * dt;
	p->vy += p->ay * dt;
	p->vz += p->az * dt;
}

