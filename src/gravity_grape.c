/**
 * @file 	gravity.c
 * @brief 	Gravity calculation using GRAPE.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	GRAPE is a special purpose hardware to accelerate
 * N-body simulations. This routine calculates self gravity using
 * a GRAPE 7 card. 
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
double gravity_minimum_mass = 1e300;

void gravity_calculate_acceleration(){
	// Initialize grape
	if (grape_open==0){
		grape_open=1;
		printf("\n************** GRAPE STARTING ***************\n");
		g5_open();
		jmemsize = g5_get_jmemsize();
		printf("************** GRAPE STARTED  ***************\n");
	}

	g5_set_range(-boxsize/2., boxsize/2., gravity_minimum_mass);
	
#ifdef INTEGRATOR_WH
	const int firstParticle = 1;
#else //INTEGRATOR_WH
	const int firstParticle = 0;
#endif //INTEGRATOR_WH
	
	// Initialize or increase memory if needed.
	int	nj = (N_active==-1)?N:N_active;		// Massive particles
	int 	ni = N;					// All particles
	nj -= firstParticle;				// Do not sum over central object for WH	
		
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
	// Copy active (j) particle mass and positions
	for(int j=0;j<nj;j++){
		mj[j] 		= particles[j+firstParticle].m;
		xj[j][0] 	= particles[j+firstParticle].x;
		xj[j][1] 	= particles[j+firstParticle].y;
		xj[j][2] 	= particles[j+firstParticle].z;
	}
	for(int i=0;i<ni;i++){
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
	}

	// Summing over all Ghost Boxes
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
	for (int gbz=-nghostz; gbz<=nghostz; gbz++){
		struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,gbz);
		// Copy passive (i) particle mass and positions
		for(int i=0;i<ni;i++){
			xi[i][0] 	= gb.shiftx+particles[i].x;
			xi[i][1] 	= gb.shifty+particles[i].y;
			xi[i][2] 	= gb.shiftz+particles[i].z;
		}
		// Summing over all particle pairs
		int nj_cur = 0;
		while (nj_cur<nj){	
			int nj_tmp = nj-nj_cur>jmemsize?jmemsize:nj-nj_cur;
			g5_set_jp(0, nj_tmp, &(mj[nj_cur]), &(xj[nj_cur])); 
			g5_set_n(nj_tmp); 
			g5_set_eps_to_all(softening);
			g5_calculate_force_on_x(xi, ai, pi, ni); 
		
			// Updateing acceleration	
			for(int i=0;i<ni;i++){
				particles[i].ax += G*ai[i][0];
				particles[i].ay += G*ai[i][1];
				particles[i].az += G*ai[i][2];
			}
			nj_cur+=nj_tmp;
		}
	}
	}
	}
}

// Try to close GRAPE
void gravity_finish(){
	if (grape_open==1){
		grape_open=0;
		g5_close();
	}
}
