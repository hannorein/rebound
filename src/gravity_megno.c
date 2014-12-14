/**
 * @file 	gravity.c
 * @brief 	Direct gravity calculation, O(N^2).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	This is the crudest implementation of an N-body code
 * which sums up every pair of particles. It is only useful very small 
 * particle numbers (N<~100) as it scales as O(N^2). Note that the MPI
 * implementation is not well tested and only works for very specific
 * problems. This should be resolved in the future. 
 *
 * 
 * @section LICENSE
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
#include "main.h"
#include "tools.h"
#include "boundaries.h"
#include "communication_mpi.h"

#ifndef INTEGRATOR_IAS15
#error GRAVITY_MEGNO requires INTERGRATOR_IAS15
#endif

extern double dt_last_success;

void gravity_calculate_acceleration(){
#pragma omp parallel for schedule(guided)
	for (int i=0; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
	// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
	for (int i=0; i<N/2; i++){
		double csx = 0;
		double csy = 0;
		double csz = 0;
	for (int j=0; j<N/2; j++){
		if (i==j) continue;
		double dx = particles[i].x - particles[j].x;
		double dy = particles[i].y - particles[j].y;
		double dz = particles[i].z - particles[j].z;
		double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
		double prefact = -G/(r*r*r)*particles[j].m;
		
		double ax = particles[i].ax;
		csx  +=	prefact*dx; 
		particles[i].ax    = ax + csx;
		csx  += ax - particles[i].ax; 
		
		double ay = particles[i].ay;
		csy  +=	prefact*dy; 
		particles[i].ay    = ay + csy;
		csy  += ay - particles[i].ay; 
		
		double az = particles[i].az;
		csz  +=	prefact*dz; 
		particles[i].az    = az + csz;
		csz  += az - particles[i].az; 
	}
	}
	/// MEGNO
#pragma omp parallel for schedule(guided)
	for (int i=N/2; i<N; i++){
	for (int j=N/2; j<N; j++){
		if (i==j) continue;
		const double dx = particles[i-N/2].x - particles[j-N/2].x;
		const double dy = particles[i-N/2].y - particles[j-N/2].y;
		const double dz = particles[i-N/2].z - particles[j-N/2].z;
		const double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
		const double r3inv = 1./(r*r*r);
		const double r5inv = 3./(r*r*r*r*r);
		const double ddx = particles[i].x - particles[j].x;
		const double ddy = particles[i].y - particles[j].y;
		const double ddz = particles[i].z - particles[j].z;
		const double Gm = G * particles[j].m;
		
		// Variational equations
		particles[i].ax += Gm * (
			+ ddx * ( dx*dx*r5inv - r3inv )
			+ ddy * ( dx*dy*r5inv )
			+ ddz * ( dx*dz*r5inv )
			);

		particles[i].ay += Gm * (
			+ ddx * ( dy*dx*r5inv )
			+ ddy * ( dy*dy*r5inv - r3inv )
			+ ddz * ( dy*dz*r5inv )
			);

		particles[i].az += Gm * (
			+ ddx * ( dz*dx*r5inv )
			+ ddy * ( dz*dy*r5inv )
			+ ddz * ( dz*dz*r5inv - r3inv )
			);
		
	}
	}
}

/// MEGNO helper routines
double gravity_megno_Ys;
double gravity_megno_Yss;
void gravity_megno_init(){
	int Nreal = N;
	gravity_megno_Ys = 0.;
	gravity_megno_Yss = 0.;
        for (int i=0;i<Nreal;i++){ 
                struct particle megno;
                megno.m = particles[i].m;
                megno.x  = 1e-16*tools_normal(1.);
                megno.y  = 1e-16*tools_normal(1.);
                megno.z  = 1e-16*tools_normal(1.);
                megno.vx = 1e-16*tools_normal(1.);
                megno.vy = 1e-16*tools_normal(1.);
                megno.vz = 1e-16*tools_normal(1.);
                particles_add(megno);
        }
}

double gravity_megno_delta(){
        double delta = 0;
        for (int i=N/2;i<N;i++){
                delta += particles[i].x * particles[i].x; 
                delta += particles[i].y * particles[i].y;
                delta += particles[i].z * particles[i].z;
                delta += particles[i].vx * particles[i].vx; 
                delta += particles[i].vy * particles[i].vy;
                delta += particles[i].vz * particles[i].vz;
        }
        return sqrt(delta);
}
double gravity_megno_deltad(){
        double deltad = 0;
        for (int i=N/2;i<N;i++){
                deltad += particles[i].vx * particles[i].x; 
                deltad += particles[i].vy * particles[i].y; 
                deltad += particles[i].vz * particles[i].z; 
                deltad += particles[i].ax * particles[i].vx; 
                deltad += particles[i].ay * particles[i].vy; 
                deltad += particles[i].az * particles[i].vz; 
        }
        return deltad;
}
void gravity_megno_update(){
	if (t<=0.) return;
	double d = gravity_megno_delta();
 	gravity_megno_Ys  += dt_last_success * t * gravity_megno_deltad()/(d*d);
	double Y = gravity_megno_Ys*2./t; 
	gravity_megno_Yss += Y * dt_last_success;
}
double gravity_megno(){
	return gravity_megno_Yss/t;
}
	
