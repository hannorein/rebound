/**
 * @file 	boundaries.c
 * @brief 	Implementation of open boundary conditions. 
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The code supports different boundary conditions.
 * This file implements open boundary conditions. Every particle that leaves 
 * the box is removed from the simulation 
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
#include "integrator.h"
#include "rebound.h"
#include "boundaries.h"
#include "tree.h"

void boundaries_check(struct Rebound* const r){
	struct Particle* const particles = r->particles;
	for (int i=0;i<r->N;i++){ // run through loop backwards so we don't have to recheck same index after removing
		int removep = 0;
		if(particles[i].x>r->boxsize_x/2.){
			removep = 1;
		}
		if(particles[i].x<-r->boxsize_x/2.){
			removep = 1;
		}
		if(particles[i].y>r->boxsize_y/2.){
			removep = 1;
		}
		if(particles[i].y<-r->boxsize_y/2.){
			removep = 1;
		}
		if(particles[i].z>r->boxsize_z/2.){
			removep = 1;
		}
		if(particles[i].z<-r->boxsize_z/2.){
			removep = 1;
		}
		if (removep==1){
#ifndef TREE
			particles_remove(r, i,0); // keepSorted=0 by default in C version
			i--; // need to recheck the particle that replaced the removed one
#endif
		}
	}
}

struct Ghostbox boundaries_get_ghostbox(struct Rebound* const r, int i, int j, int k){
	struct Ghostbox gb;
	gb.shiftx = r->boxsize_x*(double)i;
	gb.shifty = r->boxsize_y*(double)j;
	gb.shiftz = r->boxsize_z*(double)k;
	gb.shiftvx = 0;
	gb.shiftvy = 0;
	gb.shiftvz = 0;
	return gb;
}

/**
 * Checks if a given particle is within the computational domain.
 * @param p Particle to be checked.
 * @return Return value is 1 if particle is inside the box and 0 otherwise.
 */
int boundaries_particle_is_in_box(struct Rebound* const r, struct Particle p){
	if(p.x>r->boxsize_x/2.){
		return 0;
	}
	if(p.x<-r->boxsize_x/2.){
		return 0;
	}
	if(p.y>r->boxsize_y/2.){
		return 0;
	}
	if(p.y<-r->boxsize_y/2.){
		return 0;
	}
	if(p.z>r->boxsize_z/2.){
		return 0;
	}
	if(p.z<-r->boxsize_z/2.){
		return 0;
	}
	return 1;
}

