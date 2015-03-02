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
#include "integrator.h"
#include "boundaries.h"
#include "communication_mpi.h"

#ifdef MPI
#warning GRAVITY_DIRECT might not work with MPI for your problem. 
#warning Make sure you know what the code is doing. Have a look at the example restricted_threebody_mpi.
#endif

void gravity_calculate_acceleration(){
#pragma omp parallel for schedule(guided)
	for (int i=0; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
	// Summing over all Ghost Boxes
	const int _N_active = ((N_active==-1)?N:N_active)- N_megno;
	const int _N_real   = N - N_megno;
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
	for (int gbz=-nghostz; gbz<=nghostz; gbz++){
		struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,gbz);
		// Summing over all particle pairs
#pragma omp parallel for schedule(guided)
		for (int i=(selected_integrator==2?1:0); i<_N_real; i++){
			double csx = 0;
			double csy = 0;
			double csz = 0;
		for (int j=(selected_integrator==2?1:0); j<_N_active; j++){
			if (selected_integrator==1 && i==1 && j==0) continue;
			if (selected_integrator==1 && j==1 && i==0) continue;
			if (i==j) continue;
			double dx = (gb.shiftx+particles[i].x) - particles[j].x;
			double dy = (gb.shifty+particles[i].y) - particles[j].y;
			double dz = (gb.shiftz+particles[i].z) - particles[j].z;
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
	}
	}
	}
#ifdef MPI
	// Distribute active particles from root to all other nodes.
	// This assures that round-off errors do not accumulate and 
	// the copies of active particles do not diverge. 
	MPI_Bcast(particles, N_active, mpi_particle, 0, MPI_COMM_WORLD); 
#endif

}

void gravity_calculate_variational_acceleration(){
#pragma omp parallel for schedule(guided)
	for (int i=N-N_megno+0; i<N; i++){
	for (int j=N-N_megno+0; j<N; j++){
		if (selected_integrator==1 && i==N-N_megno+1 && j==N-N_megno) continue;
		if (selected_integrator==1 && j==N-N_megno+1 && i==N-N_megno) continue;
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
