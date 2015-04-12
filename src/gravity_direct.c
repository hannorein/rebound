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

#ifdef MPI
#include "communication_mpi.h"
#warning GRAVITY_DIRECT might not work with MPI for your problem. 
#warning Make sure you know what the code is doing. Have a look at the example restricted_threebody_mpi.
#endif

struct cs_3d {
	double x;
	double y;
	double z;
};


static struct cs_3d* restrict cs = NULL;
static int N_cs = 0;
unsigned int gravity_ignore_10;

void gravity_calculate_acceleration(){
	const unsigned int _gravity_ignore_10 = gravity_ignore_10;
	const int _N_real   = N - N_megno;
	if (N_cs<_N_real){
		cs = realloc(cs,_N_real*sizeof(struct cs_3d));
		N_cs = _N_real;
	}
#pragma omp parallel for schedule(guided)
	for (int i=0; i<_N_real; i++){
		particles[i].ax = 0.; 
		particles[i].ay = 0.; 
		particles[i].az = 0.; 
		cs[i].x = 0.;
		cs[i].y = 0.;
		cs[i].z = 0.;
	}
	const int _N_active = ((N_active==-1)?N:N_active)- N_megno;
	const int _N_start  = (integrator==WH?1:0);
	// Summing over all massive particle pairs
#pragma omp parallel for schedule(guided)
	for (int i=_N_start; i<_N_active; i++){
	for (int j=i+1; j<_N_active; j++){
		if (_gravity_ignore_10 && j==1 && i==0 ) continue;
		const double dx = particles[i].x - particles[j].x;
		const double dy = particles[i].y - particles[j].y;
		const double dz = particles[i].z - particles[j].z;
		const double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
		const double r = sqrt(r2);
		const double prefact  = -G/(r2*r);
		const double prefacti = prefact*particles[i].m;
		const double prefactj = prefact*particles[j].m;
		
		{
		double ax = particles[i].ax;
		cs[i].x  +=	prefactj*dx; 
		particles[i].ax    = ax + cs[i].x;
		cs[i].x  += ax - particles[i].ax; 
		
		double ay = particles[i].ay;
		cs[i].y  +=	prefactj*dy; 
		particles[i].ay    = ay + cs[i].y;
		cs[i].y  += ay - particles[i].ay; 
		
		double az = particles[i].az;
		cs[i].z  +=	prefactj*dz; 
		particles[i].az    = az + cs[i].z;
		cs[i].z  += az - particles[i].az; 
		}
		
		{
		double ax = particles[j].ax;
		cs[j].x  -=	prefacti*dx; 
		particles[j].ax    = ax + cs[j].x;
		cs[j].x  += ax - particles[j].ax; 
		
		double ay = particles[j].ay;
		cs[j].y  -=	prefacti*dy; 
		particles[j].ay    = ay + cs[j].y;
		cs[j].y  += ay - particles[j].ay; 
		
		double az = particles[j].az;
		cs[j].z  -=	prefacti*dz; 
		particles[j].az    = az + cs[j].z;
		cs[j].z  += az - particles[j].az; 
		}
		
	}
	}
	// Testparticles
#pragma omp parallel for schedule(guided)
	for (int i=_N_active; i<_N_real; i++){
	for (int j=_N_start; j<_N_active; j++){
		if (_gravity_ignore_10 && ((i==1 && j==0) || (j==1 && i==0)) ) continue;
		const double dx = particles[i].x - particles[j].x;
		const double dy = particles[i].y - particles[j].y;
		const double dz = particles[i].z - particles[j].z;
		const double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
		const double r = sqrt(r2);
		const double prefact = -G/(r2*r)*particles[j].m;
		
		double ax = particles[i].ax;
		cs[i].x  +=	prefact*dx; 
		particles[i].ax    = ax + cs[i].x;
		cs[i].x  += ax - particles[i].ax; 
		
		double ay = particles[i].ay;
		cs[i].y  +=	prefact*dy; 
		particles[i].ay    = ay + cs[i].y;
		cs[i].y  += ay - particles[i].ay; 
		
		double az = particles[i].az;
		cs[i].z  +=	prefact*dz; 
		particles[i].az    = az + cs[i].z;
		cs[i].z  += az - particles[i].az; 
	}
	}
		
}

void gravity_calculate_variational_acceleration(){
	const unsigned int _gravity_ignore_10 = gravity_ignore_10;
	const int _N_real   = N - N_megno;
#pragma omp parallel for schedule(guided)
	for (int i=_N_real; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
#pragma omp parallel for schedule(guided)
	for (int i=_N_real; i<N; i++){
	for (int j=i+1; j<N; j++){
		if (_gravity_ignore_10 && ((i==_N_real+1 && j==_N_real) || (j==_N_real+1 && i==_N_real)) ) continue;
		const double dx = particles[i-N/2].x - particles[j-N/2].x;
		const double dy = particles[i-N/2].y - particles[j-N/2].y;
		const double dz = particles[i-N/2].z - particles[j-N/2].z;
		const double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
		const double r  = sqrt(r2);
		const double r3inv = 1./(r2*r);
		const double r5inv = 3.*r3inv/r2;
		const double ddx = particles[i].x - particles[j].x;
		const double ddy = particles[i].y - particles[j].y;
		const double ddz = particles[i].z - particles[j].z;
		const double Gmi = G * particles[i].m;
		const double Gmj = G * particles[j].m;
		
		// Variational equations
		const double dax =   ddx * ( dx*dx*r5inv - r3inv )
				   + ddy * ( dx*dy*r5inv )
				   + ddz * ( dx*dz*r5inv );
		const double day =   ddx * ( dy*dx*r5inv )
				   + ddy * ( dy*dy*r5inv - r3inv )
				   + ddz * ( dy*dz*r5inv );
		const double daz =   ddx * ( dz*dx*r5inv )
				   + ddy * ( dz*dy*r5inv )
				   + ddz * ( dz*dz*r5inv - r3inv );
		
		particles[i].ax += Gmj * dax;
		particles[i].ay += Gmj * day;
		particles[i].az += Gmj * daz;
		
		particles[j].ax -= Gmi * dax;
		particles[j].ay -= Gmi * day;
		particles[j].az -= Gmi * daz;
	}
	}
}
