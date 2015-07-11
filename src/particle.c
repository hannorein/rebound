/**
 * @file 	particle.c
 * @brief 	Particle structure and main particle routines.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
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
#include "tree.h"
#ifndef COLLISIONS_NONE
#include "collisions.h"
#endif // COLLISIONS_NONE
#ifdef MPI
#include "communication_mpi.h"
#endif // MPI

struct particle* 	particles = NULL;	

#ifdef BOUNDARIES_OPEN
int boundaries_particle_is_in_box(struct particle p);
int particles_warning_is_not_in_box = 0;
#endif // BOUNDARIES_OPEN
#ifdef GRAVITY_GRAPE
extern double gravity_minimum_mass;
#endif // GRAVITY_GRAPE

void particles_add_local(struct Rebound* const r, struct particle pt){
#ifdef BOUNDARIES_OPEN
	if (boundaries_particle_is_in_box(pt)==0){
		// Particle has left the box. Do not add.
		if (particles_warning_is_not_in_box == 0){
			particles_warning_is_not_in_box = 1;
			printf("\nWarning: Trying to add particle which is outside box boundaries.\n");
		}
		return;
	}
#endif // BOUNDARIES_OPEN
	while (r->Nmax<=r->N){
		r->Nmax += 128;
		particles = realloc(particles,sizeof(struct particle)*r->Nmax);
	}

	particles[r->N] = pt;

#ifdef TREE
	tree_add_particle_to_tree(r->N);
#endif // TREE
	(r->N)++;
}

void particles_add(struct Rebound* const r, struct particle pt){
	if (r->N_megno){
		printf("\nWarning: Trying to add particle after calling megno_init().\n");
	}
#ifndef COLLISIONS_NONE
	if (pt.r>=collisions_max_r){
		collisions_max2_r = collisions_max_r;
		collisions_max_r = pt.r;
	}else{
		if (pt.r>=collisions_max2_r){
			collisions_max2_r = pt.r;
		}
	}
#endif 	// COLLISIONS_NONE
#ifdef GRAVITY_GRAPE
	if (pt.m<gravity_minimum_mass){
		gravity_minimum_mass = pt.m;
	}
#endif // GRAVITY_GRAPE
#ifdef MPI
	int rootbox = particles_get_rootbox_for_particle(r, pt);
	int root_n_per_node = root_n/mpi_num;
	int proc_id = rootbox/root_n_per_node;
	if (proc_id != mpi_id && r->N >= r->N_active){
		// Add particle to array and send them to proc_id later. 
		communication_mpi_add_particle_to_send_queue(pt,proc_id);
		return;
	}
#endif // MPI
	// Add particle to local partical array.
	particles_add_local(r, pt);
}

void particles_add_fixed(struct Rebound* const r, struct particle pt,int pos){
	// Only works for non-MPI simulations or when the particles does not move to another node.
#ifdef BOUNDARIES_OPEN
	if (boundaries_particle_is_in_box(pt)==0){
		// Particle has left the box. Do not add.
		return;
	}
#endif // BOUNDARIES_OPEN
	particles[pos] = pt; 
#ifdef TREE
	tree_add_particle_to_tree(pos);
#endif // TREE
}

int particles_get_rootbox_for_particle(const struct Rebound* const r, struct particle pt){
	if (r->boxsize==-1) return 0;
	int i = ((int)floor((pt.x + r->boxsize_x/2.)/r->boxsize)+r->root_nx)%r->root_nx;
	int j = ((int)floor((pt.y + r->boxsize_y/2.)/r->boxsize)+r->root_ny)%r->root_ny;
	int k = ((int)floor((pt.z + r->boxsize_z/2.)/r->boxsize)+r->root_nz)%r->root_nz;
	int index = (k*r->root_ny+j)*r->root_nx+i;
	return index;
}

void particles_remove_all(struct Rebound* const r){
	r->N 		= 0;
	r->Nmax 	= 0;
	r->N_active 	= -1;
	r->N_megno 	= 0;
	free(particles);
	particles 	= NULL;
}

int particles_remove(struct Rebound* const r, int index, int keepSorted){
	if (r->N==1){
		fprintf(stderr, "Last particle removed.\n");
		return 1;
	}
	if (index >= r->N){
		fprintf(stderr, "\nIndex %d passed to particles_remove was out of range (N=%d).  Did not remove particle.\n", index, r->N);
		return 0;
	}
	if (r->N_megno){
		fprintf(stderr, "\nRemoving particles not supported when calculating MEGNO.  Did not remove particle.\n");
		return 0;
	}
	(r->N)--;
	if(keepSorted){
		for(int j=index; j<r->N; j++){
			particles[j] = particles[j+1];
		}
	}
	else{
		particles[index] = particles[r->N];
	}

	return 1;
}

#ifdef PARTICLEIDS
int particles_remove_ID(struct Rebound* const r, int ID, int keepSorted){
	int success = 0;
	for(int i=0;i<r->N;i++){
		if(particles[i].ID == ID){
			success = particles_remove(i, keepSorted);
			break;
		}
	}

	if(!success){
		fprintf(stderr, "\nIndex passed to particles_remove_ID (ID = %d) not found in particles array.  Did not remove particle.\n", ID);
	}
	return success;
}
#endif
