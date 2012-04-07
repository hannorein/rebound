/**
 * @file 	collisions.c
 * @brief 	Collision search using an oct-tree, O(N log(N)).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details 	The routines in this file implement a nearest neighbour 
 * search using the oct-tree defined in file tree.h. It can be run with
 * MPI using a distributed tree. In that case the locally essential tree 
 * is shared with every other node. The method scales as O(N log(N)) for 
 * large particles. For small particle numbers, a direct collision
 * search might be faster, as it avoids having the overhead of a 
 * complicated data structure. 
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
#include "tree.h"
#include "boundaries.h"
#include "communication_mpi.h"

struct 	collision* collisions 	= NULL;		/**< Array of all collisions. */
int 	collisions_NMAX 	= 0;		/**< Size allocated for collisions.*/
int 	collisions_N 		= 0;		/**< Number of elements in collisions. */

double 	collisions_max_r	= 0;
double 	collisions_max2_r	= 0;

/**
 * Find the nearest neighbour in a cell or its daughters.
 * The function only returns a positive result if the particles
 * are overlapping. Thus, the name nearest neighbour is not
 * exactly true.
 * @param gb (Shifted) position and velocity of the particle.
 * @param ri Index of the root box currently being searched in.
 * @param p1_r Radius of the particle (this is not in gb).
 * @param nearest_r2 Pointer to the nearest neighbour found so far.
 * @param collision_nearest Pointer to the nearest collision found so far.
 * @param c Pointer to the cell currently being searched in.
 */
void tree_get_nearest_neighbour_in_cell(struct ghostbox gb, struct ghostbox gbunmod, int ri, double p1_r,  double* nearest_r2, struct collision* collision_nearest, struct cell* c);

void collisions_search(){
	// Update and simplify tree. 
	// Prepare particles for distribution to other nodes. 
	tree_update();          

#ifdef MPI
	// Distribute particles and add newly received particles to tree.
	communication_mpi_distribute_particles();
	
	// Prepare essential tree (and particles close to the boundary needed for collisions) for distribution to other nodes.
	tree_prepare_essential_tree_for_collisions();

	// Transfer essential tree and particles needed for collisions.
	communication_mpi_distribute_essential_tree_for_collisions();
#endif // MPI

	// Loop over ghost boxes, but only the inner most ring.
	int nghostxcol = (nghostx>1?1:nghostx);
	int nghostycol = (nghosty>1?1:nghosty);
	int nghostzcol = (nghostz>1?1:nghostz);
	// Loop over all particles
#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		struct particle p1 = particles[i];
		struct  collision collision_nearest;
		collision_nearest.p1 = i;
		collision_nearest.p2 = -1;
		double p1_r = p1.r;
		double nearest_r2 = boxsize_max*boxsize_max/4.;
		// Loop over ghost boxes.
		for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
		for (int gby=-nghostycol; gby<=nghostycol; gby++){
		for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
			// Calculated shifted position (for speedup). 
			struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,gbz);
			struct ghostbox gbunmod = gb;
			gb.shiftx += p1.x; 
			gb.shifty += p1.y; 
			gb.shiftz += p1.z; 
			gb.shiftvx += p1.vx; 
			gb.shiftvy += p1.vy; 
			gb.shiftvz += p1.vz; 
			// Loop over all root boxes.
			for (int ri=0;ri<root_n;ri++){
				struct cell* rootcell = tree_root[ri];
				if (rootcell!=NULL){
					tree_get_nearest_neighbour_in_cell(gb, gbunmod,ri,p1_r,&nearest_r2,&collision_nearest,rootcell);
				}
			}
		}
		}
		}
		// Continue if no collision was found
		if (collision_nearest.p2==-1) continue;
	}
}

/**
 * This is really nice with all those parametes.
 * @TODO Cleanup
 */ 
void tree_get_nearest_neighbour_in_cell(struct ghostbox gb, struct ghostbox gbunmod, int ri, double p1_r, double* nearest_r2, struct collision* collision_nearest, struct cell* c){
	if (c->pt>=0){ 	
		// c is a leaf node
		int condition 	= 1;
#ifdef MPI
		int isloc	= 1 ;
		isloc = communication_mpi_rootbox_is_local(ri);
		if (isloc==1){
#endif // MPI
			/**
			 * If this is a local cell, make sure particle is not colliding with itself.
			 * If this is a remote cell, the particle number might be the same, even for 
			 * different particles. 
			 * @TODO This can probably be written in a cleaner way.
			 */
			condition = (c->pt != collision_nearest->p1);
#ifdef MPI
		}
#endif // MPI
		if (condition){
			struct particle p2;
#ifdef MPI
			if (isloc==1){
#endif // MPI
				p2 = particles[c->pt];
#ifdef MPI
			}else{
				int root_n_per_node = root_n/mpi_num;
				int proc_id = ri/root_n_per_node;
				p2 = particles_recv[proc_id][c->pt];
			}
#endif // MPI

			double dx = gb.shiftx - p2.x;
			double dy = gb.shifty - p2.y;
			double dz = gb.shiftz - p2.z;
			double r2 = dx*dx+dy*dy+dz*dz;
			// A closer neighbour has already been found 
			//if (r2 > *nearest_r2) return;
			double rp = p1_r+p2.r;
			// Particles are not overlapping 
			if (r2 > rp*rp) return;
			double dvx = gb.shiftvx - p2.vx;
			double dvy = gb.shiftvy - p2.vy;
			double dvz = gb.shiftvz - p2.vz;
			// Particles are not approaching each other
			if (dvx*dx + dvy*dy + dvz*dz >0) return;
			// Found a new nearest neighbour. Save it for later.
			*nearest_r2 = r2;
			collision_nearest->ri = ri;
			collision_nearest->p2 = c->pt;
			collision_nearest->gb = gbunmod;
			// Save collision in collisions array.
#pragma omp critical
			{
				if (collisions_NMAX<=collisions_N){
					collisions_NMAX += 32;
					collisions = realloc(collisions,sizeof(struct collision)*collisions_NMAX);
				}
				collisions[collisions_N] = *collision_nearest;
				collisions_N++;
			}
		}
	}else{		
		// c is not a leaf node
		double dx = gb.shiftx - c->x;
		double dy = gb.shifty - c->y;
		double dz = gb.shiftz - c->z;
		double r2 = dx*dx + dy*dy + dz*dz;
		double rp  = p1_r + collisions_max2_r + 0.86602540378443*c->w;
		// Check if we need to decent into daughter cells
		if (r2 < rp*rp ){
			for (int o=0;o<8;o++){
				struct cell* d = c->oct[o];
				if (d!=NULL){
					tree_get_nearest_neighbour_in_cell(gb,gbunmod,ri,p1_r,nearest_r2,collision_nearest,d);
				}
			}
		}
	}
}

void collisions_resolve(){
	// randomize
	for (int i=0;i<collisions_N;i++){
		int new = rand()%collisions_N;
		struct collision c1 = collisions[i];
		collisions[i] = collisions[new];
		collisions[new] = c1;
	}
	// Loop over all collisions previously found in collisions_search().
	for (int i=0;i<collisions_N;i++){
		// Resolve collision
		collision_resolve(collisions[i]);
	}
	// Mark all collisions as resolved.
	collisions_N=0;
}
