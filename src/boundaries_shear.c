/**
 * @file 	boundaries.c
 * @brief 	Implementation of shear periodic boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * This file implements shear periodic boundary conditions which are often 
 * used for studying astrophysical discs and rings. 
 * The shear is linear in the x direction (radial). The azimuthal direction
 * is y and the vertical direction is z. The  orbtial (epicyclic) frequency 
 * is set by the constant OMEGA (default: 1), which can be set in the function
 * problem_init(). It is also possible to set a different vertical epicyclic 
 * frequency with the variable OMEGAZ. For simplicity, the boundary condition
 * is periodic in the z direction. 
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
#include "boundaries.h"
#include "rebound.h"
#include "tree.h"
#include "communication_mpi.h"

int nghostx = 1;
int nghosty = 1;
int nghostz = 0;	/**< The boundary condition is periodic in z, but usually we don't need any ghostboxed as the disc is stratified */

void boundaries_check(struct Rebound* const r){
	// The offset of ghostcell is time dependent.
	const double OMEGA = r->ri_sei.OMEGA;
	double offsetp1 = -fmod(-1.5*OMEGA*r->boxsize_x*r->t+r->boxsize_y/2.,r->boxsize_y)-r->boxsize_y/2.; 
	double offsetm1 = -fmod( 1.5*OMEGA*r->boxsize_x*r->t-r->boxsize_y/2.,r->boxsize_y)+r->boxsize_y/2.; 
	struct Particle* const particles = r->particles;
#pragma omp parallel for schedule(guided)
	for (int i=0;i<r->N;i++){
		// Radial
		while(particles[i].x>r->boxsize_x/2.){
			particles[i].x -= r->boxsize_x;
			particles[i].y += offsetp1;
			particles[i].vy += 3./2.*OMEGA*r->boxsize_x;
		}
		while(particles[i].x<-r->boxsize_x/2.){
			particles[i].x += r->boxsize_x;
			particles[i].y += offsetm1;
			particles[i].vy -= 3./2.*OMEGA*r->boxsize_x;
		}
		// Azimuthal
		while(particles[i].y>r->boxsize_y/2.){
			particles[i].y -= r->boxsize_y;
		}
		while(particles[i].y<-r->boxsize_y/2.){
			particles[i].y += r->boxsize_y;
		}
		// Vertical (there should be no boundary, but periodic makes life easier)
		while(particles[i].z>r->boxsize_z/2.){
			particles[i].z -= r->boxsize_z;
		}
		while(particles[i].z<-r->boxsize_z/2.){
			particles[i].z += r->boxsize_z;
		}
	}
}

struct Ghostbox boundaries_get_ghostbox(struct Rebound* const r, int i, int j, int k){
	const double OMEGA = r->ri_sei.OMEGA;
	struct Ghostbox gb;
	// Ghostboxes habe a finite velocity.
	gb.shiftvx = 0;
	gb.shiftvy = -1.5*(double)i*OMEGA*r->boxsize_x;
	gb.shiftvz = 0;
	// The shift in the y direction is time dependent. 
	double shift;
	if (i==0){
		shift = -fmod(gb.shiftvy*r->t,r->boxsize_y); 
	}else{
		if (i>0){
			shift = -fmod(gb.shiftvy*r->t-r->boxsize_y/2.,r->boxsize_y)-r->boxsize_y/2.; 
		}else{
			shift = -fmod(gb.shiftvy*r->t+r->boxsize_y/2.,r->boxsize_y)+r->boxsize_y/2.; 
		}	
	}
	gb.shiftx = r->boxsize_x*(double)i;
	gb.shifty = r->boxsize_y*(double)j-shift;
	gb.shiftz = r->boxsize_z*(double)k;
	return gb;
}


