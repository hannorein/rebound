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
#include "boundaries.h"
#include "communication_mpi.h"

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
		const double r5inv = 1./(r*r*r*r*r);
		const double ddx = particles[i].x - particles[j].x;
		const double ddy = particles[i].y - particles[j].y;
		const double ddz = particles[i].z - particles[j].z;
		const double Gm = G * particles[j].m;
		
		// Variational equations
		particles[i].ax += Gm * (
			+ ddx * ( + 3.*dx*dx*r5inv - 1.*r3inv)
			+ ddy * ( + 3.*dx*dy*r5inv)
			+ ddz * ( + 3.*dx*dz*r5inv)
			);

		particles[i].ay += Gm * (
			+ ddx * ( + 3.*dy*dx*r5inv)
			+ ddy * ( + 3.*dy*dy*r5inv - 1.*r3inv)
			+ ddz * ( + 3.*dy*dz*r5inv)
			);

		particles[i].az += Gm * (
			+ ddx * ( + 3.*dz*dx*r5inv)
			+ ddy * ( + 3.*dz*dy*r5inv)
			+ ddz * ( + 3.*dz*dz*r5inv - 1.*r3inv)
			);
		
	}
	}
}
