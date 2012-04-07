/**
 * @file 	collisions.c
 * @brief 	Direct collision search, O(N^2).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	This is the crudest implementation of a collision search
 * and checks every pair of particles. It is only useful very small 
 * particle numbers (N<~100) as it scales as O(N^2). 
 *
 * A collision is defined as an overlap between two particles. This
 * is only an approximation and works only if the timestep is small
 * enough. More precisely, dt << v / Rp, where v is the typical velocity
 * and Rp the radius of a particle. Furthermore, particles must be 
 * approaching each other at the time when they overlap. 
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
#include "collisions.h"
#include "collision_resolve.h"
#include "main.h"
#include "boundaries.h"

#ifdef MPI
#error COLLISIONS_DIRECT not yet compatible with MPI
#endif

struct	collision* collisions 	= NULL;		/**< Array of all collisions. */
int 	collisions_NMAX 	= 0;		/**< Size allocated for collisions.*/
int 	collisions_N 		= 0;		/**< Number of elements in collisions. */

double 	collisions_max_r	= 0;
double 	collisions_max2_r	= 0;

void collisions_search(){
	// Loop over ghost boxes, but only the inner most ring.
	int nghostxcol = (nghostx>1?1:nghostx);
	int nghostycol = (nghosty>1?1:nghosty);
	int nghostzcol = (nghostz>1?1:nghostz);
	for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
	for (int gby=-nghostycol; gby<=nghostycol; gby++){
	for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
		// Loop over all particles
		for (int i=0;i<N;i++){
			struct particle p1 = particles[i];
			struct ghostbox gborig = boundaries_get_ghostbox(gbx,gby,gbz);
			struct ghostbox gb = gborig;
			// Precalculate shifted position 
			gb.shiftx += p1.x;
			gb.shifty += p1.y;
			gb.shiftz += p1.z;
			gb.shiftvx += p1.vx;
			gb.shiftvy += p1.vy;
			gb.shiftvz += p1.vz;
			// Loop over all particles again
			for (int j=0;j<N;j++){
				// Do not collide particle with itself.
				if (i==j) continue;
				struct particle p2 = particles[j];
				double dx = gb.shiftx - p2.x; 
				double dy = gb.shifty - p2.y; 
				double dz = gb.shiftz - p2.z; 
				double sr = p1.r + p2.r; 
				double r2 = dx*dx+dy*dy+dz*dz;
				// Check if particles are overlapping 
				if (r2>sr*sr) continue;	
				double dvx = gb.shiftvx - p2.vx; 
				double dvy = gb.shiftvy - p2.vy; 
				double dvz = gb.shiftvz - p2.vz; 
				// Check if particles are approaching each other
				if (dvx*dx + dvy*dy + dvz*dz >0) continue; 
				// Add particles to collision array.
				if (collisions_NMAX<=collisions_N){
					// Allocate memory if there is no space in array.
					// Doing it in chunks of 32 to avoid having to do it too often.
					collisions_NMAX += 32;
					collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
				}
				collisions[collisions_N].p1 = i;
				collisions[collisions_N].p2 = j;
				collisions[collisions_N].gb = gborig;
				collisions_N++;
			}
		}
	}
	}
	}
}

void collisions_resolve(){
	// Loop over all collisions previously found in collisions_search().
	for (int i=0;i<collisions_N;i++){
		// Resolve collision
		collision_resolve(collisions[i]);
	}
	// Mark all collisions as resolved.
	collisions_N=0;
}
