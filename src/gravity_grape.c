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
#include "g5nbutil.h"
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "communication_mpi.h"

int _nj_MAX 	= 0;
int _ni_MAX 	= 0;
double* mj	= NULL;
double* pi	= NULL;
double (*xj)[3]	= NULL;
double (*xi)[3]	= NULL;	
double (*ai)[3]	= NULL;


int grape_open	= 0;
int jmemsize 	= 0;

void gravity_calculate_acceleration(){
	// Initialize grape
	if (grape_open==0){
		grape_open=1;
		printf("\n************** GRAPE START ***************\n");
		g5_open();
		jmemsize = g5_get_jmemsize();
		printf("************** GRAPE START ***************\n");
	}

	g5_set_range(-boxsize/2., boxsize/2., 1e-6);
	
	// Initialize or increase memory if needed.
	int	nj = (N_active==-1)?N:N_active;		// Massive particles
	int 	ni = N;					// All particles
	if (_nj_MAX<nj){
		_nj_MAX = nj;
		mj 	= realloc(mj,sizeof(double)*nj);
		xj	= realloc(xj,sizeof(double)*nj*3);
	}
	if (_ni_MAX<ni){
		_ni_MAX = ni;
		xi	= realloc(xi,sizeof(double)*ni*3);	// position
		ai	= realloc(ai,sizeof(double)*ni*3);	// acceleration
		pi 	= realloc(pi,sizeof(double)*ni);	// potential
	}
	// Copy particle mass and positions
	for(int i=0;i<nj;i++){
		mj[i] 		= particles[i].m;
		xj[i][0] 	= particles[i].x;
		xj[i][1] 	= particles[i].y;
		xj[i][2] 	= particles[i].z;
	}
	for(int i=0;i<ni;i++){
		xi[i][0] 	= particles[i].x;
		xi[i][1] 	= particles[i].y;
		xi[i][2] 	= particles[i].z;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
	}

	int nj_cur = 0;
	while (nj_cur<nj){	
		int nj_tmp = nj-nj_cur>jmemsize?jmemsize:nj-nj_cur;
		g5_set_jp(0, nj_tmp, &(mj[nj_cur]), &(xj[nj_cur])); 
		g5_set_n(nj_tmp); 
		g5_set_eps_to_all(softening);
		g5_calculate_force_on_x(xi, ai, pi, ni); 
		
		for(int i=0;i<ni;i++){
			particles[i].ax += ai[i][0];
			particles[i].ay += ai[i][1];
			particles[i].az += ai[i][2];
		}
		nj_cur+=nj_tmp;
	}
}

void gravity_finish(){
	if (grape_open==1){
		grape_open=0;
		g5_close();
	}
}
