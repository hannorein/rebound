/**
 * @file 	communication_mpi.c
 * @brief	Handles communication between nodes using MPI.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 *
 * @details These routines handle the communication between
 * different nodes via the Message Passing Interface (MPI). 
 * There are two different types of communications implemented 
 * at the moment:
 * - Distributing particles to the correct node.
 * - Creating, and distributing the essential tree to allow
 *   other nodes walk remote trees. Note that the opening 
 *   criteria is different for gravity and collision 
 *   tree walks.
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
#ifdef MPI
#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "rebound.h"
#include "tree.h"
#include "boundary.h"
#include "communication_mpi.h"

void reb_communication_mpi_init(struct reb_simulation* const r, int argc, char** argv){
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized){
	    MPI_Init(&argc,&argv);
    }
	MPI_Comm_size(MPI_COMM_WORLD,&(r->mpi_num));
	MPI_Comm_rank(MPI_COMM_WORLD,&(r->mpi_id));
	
	// Prepare send/recv buffers for particles
	r->particles_send   	= calloc(r->mpi_num,sizeof(struct reb_particle*));
	r->N_particles_send 	= calloc(r->mpi_num,sizeof(int));
	r->N_particles_send_max 	= calloc(r->mpi_num,sizeof(int));
	r->particles_recv   	= calloc(r->mpi_num,sizeof(struct reb_particle*));
	r->N_particles_recv 	= calloc(r->mpi_num,sizeof(int));
	r->N_particles_recv_max 	= calloc(r->mpi_num,sizeof(int));

	// Prepare send/recv buffers for essential tree
	r->tree_essential_send   	= calloc(r->mpi_num,sizeof(struct reb_treecell*));
	r->N_tree_essential_send 	= calloc(r->mpi_num,sizeof(int));
	r->N_tree_essential_send_max = calloc(r->mpi_num,sizeof(int));
	r->tree_essential_recv   	= calloc(r->mpi_num,sizeof(struct reb_treecell*));
	r->N_tree_essential_recv 	= calloc(r->mpi_num,sizeof(int));
	r->N_tree_essential_recv_max = calloc(r->mpi_num,sizeof(int));
}

int reb_communication_mpi_rootbox_is_local(struct reb_simulation* const r, int i){
	int N_root_per_node = r->N_root/r->mpi_num;
	int proc_id = i/N_root_per_node;
	if (proc_id != r->mpi_id){
		return 0;
	}else{
		return 1;
	}
}


void reb_communication_mpi_distribute_particles(struct reb_simulation* const r){
	// Distribute the number of particles to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->N_particles_send, 1, MPI_INT, &(r->N_particles_recv[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->N_particles_recv_max[i]<r->N_particles_recv[i]){
			r->N_particles_recv_max[i] += 32;
			r->particles_recv[i] = realloc(r->particles_recv[i],sizeof(struct reb_particle)*r->N_particles_recv_max[i]);
		}
	}

	// Exchange particles via MPI.
	// Using non-blocking receive call.
	MPI_Request request[r->mpi_num];
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_particles_recv[i]==0) continue;
		MPI_Irecv(r->particles_recv[i], sizeof(struct reb_particle)*r->N_particles_recv[i], MPI_CHAR, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_particles_send[i]==0) continue;
		MPI_Send(r->particles_send[i], sizeof(struct reb_particle)* r->N_particles_send[i], MPI_CHAR, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_particles_recv[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add particles to local tree
	for (int i=0;i<r->mpi_num;i++){
		for (int j=0;j<r->N_particles_recv[i];j++){
			reb_simulation_add(r,r->particles_recv[i][j]);
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->N_particles_send[i] = 0;
		r->N_particles_recv[i] = 0;
	}
}

void reb_communication_mpi_add_particle_to_send_queue(struct reb_simulation* const r, struct reb_particle pt, int proc_id){
	int send_N = r->N_particles_send[proc_id];
	while (r->N_particles_send_max[proc_id] <= send_N){
		r->N_particles_send_max[proc_id] += 128;
		r->particles_send[proc_id] = realloc(r->particles_send[proc_id],sizeof(struct reb_particle)*r->N_particles_send_max[proc_id]);
	}
	r->particles_send[proc_id][send_N] = pt;
	r->N_particles_send[proc_id]++;
}


/** 
 * This is the data structure for an axis aligned bounding box.
 */
struct reb_aabb{
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
};

struct reb_aabb communication_boundingbox_for_root(struct reb_simulation* const r, int index){
	int i = index%r->N_root_x;
	int j = ((index-i)/r->N_root_x)%r->N_root_y;
	int k = ((index-i)/r->N_root_x-j)/r->N_root_y;
	struct reb_aabb boundingbox;
	boundingbox.xmin = -r->boxsize.x/2.+r->root_size*(double)i;
	boundingbox.ymin = -r->boxsize.y/2.+r->root_size*(double)j;
	boundingbox.zmin = -r->boxsize.z/2.+r->root_size*(double)k;
	boundingbox.xmax = -r->boxsize.x/2.+r->root_size*(double)(i+1);
	boundingbox.ymax = -r->boxsize.y/2.+r->root_size*(double)(j+1);
	boundingbox.zmax = -r->boxsize.z/2.+r->root_size*(double)(k+1);
	return boundingbox;
}

struct reb_aabb reb_communication_boundingbox_for_proc(struct reb_simulation* const r, int proc_id){
	int N_root_per_node = r->N_root/r->mpi_num;
	int root_start = proc_id*N_root_per_node;
	int root_stop  = (proc_id+1)*N_root_per_node;
	struct reb_aabb boundingbox = communication_boundingbox_for_root(r, root_start);
	for (int i=root_start+1;i<root_stop;i++){
		struct reb_aabb boundingbox2 = communication_boundingbox_for_root(r,i);
		if (boundingbox.xmin > boundingbox2.xmin) boundingbox.xmin = boundingbox2.xmin;
		if (boundingbox.ymin > boundingbox2.ymin) boundingbox.ymin = boundingbox2.ymin;
		if (boundingbox.zmin > boundingbox2.zmin) boundingbox.zmin = boundingbox2.zmin;
		if (boundingbox.xmax < boundingbox2.xmax) boundingbox.xmax = boundingbox2.xmax;
		if (boundingbox.ymax < boundingbox2.ymax) boundingbox.ymax = boundingbox2.ymax;
		if (boundingbox.zmax < boundingbox2.zmax) boundingbox.zmax = boundingbox2.zmax;
	}
	return boundingbox;
}

double reb_communication_distance2_of_aabb_to_cell(struct reb_aabb bb, struct reb_treecell* node){
	double distancex = fabs(node->x - (bb.xmin+bb.xmax)/2.)  -  (node->w + bb.xmax-bb.xmin)/2.;
	double distancey = fabs(node->y - (bb.ymin+bb.ymax)/2.)  -  (node->w + bb.ymax-bb.ymin)/2.;
	double distancez = fabs(node->z - (bb.zmin+bb.zmax)/2.)  -  (node->w + bb.zmax-bb.zmin)/2.;
	if (distancex<0) distancex =0;
	if (distancey<0) distancey =0;
	if (distancez<0) distancez =0;
	return distancex*distancex + distancey*distancey + distancez*distancez;
}

double reb_communication_distance2_of_proc_to_node(struct reb_simulation* const r, int proc_id, struct reb_treecell* node){
	int N_ghost_xcol = (r->N_ghost_x>0?1:0);
	int N_ghost_ycol = (r->N_ghost_y>0?1:0);
	int N_ghost_zcol = (r->N_ghost_z>0?1:0);
	double distance2 = r->root_size*(double)r->N_root; // A conservative estimate for the minimum distance.
	distance2 *= distance2;
	for (int gbx=-N_ghost_xcol; gbx<=N_ghost_xcol; gbx++){
	for (int gby=-N_ghost_ycol; gby<=N_ghost_ycol; gby++){
	for (int gbz=-N_ghost_zcol; gbz<=N_ghost_zcol; gbz++){
		struct reb_vec6d gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
		struct reb_aabb boundingbox = reb_communication_boundingbox_for_proc(r, proc_id);
		boundingbox.xmin+=gb.x;
		boundingbox.xmax+=gb.x;
		boundingbox.ymin+=gb.y;
		boundingbox.ymax+=gb.y;
		boundingbox.zmin+=gb.z;
		boundingbox.zmax+=gb.z;
		// calculate distance
		double distance2new = reb_communication_distance2_of_aabb_to_cell(boundingbox,node);
		if (distance2 > distance2new) distance2 = distance2new;
	}
	}
	}
	return distance2;
}

void reb_communication_mpi_prepare_essential_cell_for_collisions_for_proc(struct reb_simulation* const r, struct reb_treecell* node, int proc){
	// Add essential cell to tree_essential_send
	if (r->N_tree_essential_send[proc]>=r->N_tree_essential_send_max[proc]){
		r->N_tree_essential_send_max[proc] += 32;
		r->tree_essential_send[proc] = realloc(r->tree_essential_send[proc],sizeof(struct reb_treecell)*r->N_tree_essential_send_max[proc]);
	}
	// Copy node to send buffer
	r->tree_essential_send[proc][r->N_tree_essential_send[proc]] = (*node);
	r->N_tree_essential_send[proc]++;
	if (node->pt>=0){ // Is leaf
		// Also transmit particle (Here could be another check if the particle actually overlaps with the other box)
		if (r->N_particles_send[proc]>=r->N_particles_send_max[proc]){
			r->N_particles_send_max[proc] += 32;
			r->particles_send[proc] = realloc(r->particles_send[proc],sizeof(struct reb_treecell)*r->N_particles_send_max[proc]);
		}
		// Copy particle to send buffer
		r->particles_send[proc][r->N_particles_send[proc]] = r->particles[node->pt];
		// Update reference from cell to particle
		r->tree_essential_send[proc][r->N_tree_essential_send[proc]-1].pt = r->N_particles_send[proc];
		r->N_particles_send[proc]++;
	}else{		// Not a leaf. Check if we need to transfer daughters.
		double distance2 = reb_communication_distance2_of_proc_to_node(r, proc,node);
		double rp  = 2.*r->max_radius0 + 0.86602540378443*node->w;
		if (distance2 < rp*rp ){
			for (int o=0;o<8;o++){
				struct reb_treecell* d = node->oct[o];
				if (d==NULL) continue;
				reb_communication_mpi_prepare_essential_cell_for_collisions_for_proc(r, d,proc);
			}
		}
	}
}
void reb_communication_mpi_prepare_essential_tree_for_collisions(struct reb_simulation* const r, struct reb_treecell* root){
	if (root==NULL) return;
	// Find out which cells are needed by every other node
	for (int i=0; i<r->mpi_num; i++){
		if (i==r->mpi_id) continue;
		reb_communication_mpi_prepare_essential_cell_for_collisions_for_proc(r, root,i);	
	}
}


void reb_communication_mpi_prepare_essential_cell_for_gravity_for_proc(struct reb_simulation* const r, struct reb_treecell* node, int proc){
	// Add essential cell to tree_essential_send
	if (r->N_tree_essential_send[proc]>=r->N_tree_essential_send_max[proc]){
		r->N_tree_essential_send_max[proc] += 32;
		r->tree_essential_send[proc] = realloc(r->tree_essential_send[proc],sizeof(struct reb_treecell)*r->N_tree_essential_send_max[proc]);
	}
	// Copy node to send buffer
	r->tree_essential_send[proc][r->N_tree_essential_send[proc]] = (*node);
	r->N_tree_essential_send[proc]++;
	if (node->pt<0){		// Not a leaf. Check if we need to transfer daughters.
		double width = node->w;
		double distance2 = reb_communication_distance2_of_proc_to_node(r, proc,node);
		if ( width*width > r->opening_angle2*distance2) {
			for (int o=0;o<8;o++){
				struct reb_treecell* d = node->oct[o];
				if (d!=NULL){
					reb_communication_mpi_prepare_essential_cell_for_gravity_for_proc(r, d,proc);
				}
			}
		}
	}
}

void reb_communication_mpi_prepare_essential_tree_for_gravity(struct reb_simulation* const r,struct reb_treecell* root){
	if (root==NULL) return;
	// Find out which cells are needed by every other node
	for (int i=0; i<r->mpi_num; i++){
		if (i==r->mpi_id) continue;
		reb_communication_mpi_prepare_essential_cell_for_gravity_for_proc(r, root,i);	
	}
}

void reb_communication_mpi_distribute_essential_tree_for_gravity(struct reb_simulation* const r){
	///////////////////////////////////////////////////////////////
	// Distribute essential tree needed for gravity and collisions
	///////////////////////////////////////////////////////////////
	
	// Distribute the number of cells to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->N_tree_essential_send, 1, MPI_INT, &(r->N_tree_essential_recv[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->N_tree_essential_recv_max[i]<r->N_tree_essential_recv[i]){
			r->N_tree_essential_recv_max[i] += 32;
			r->tree_essential_recv[i] = realloc(r->tree_essential_recv[i],sizeof(struct reb_treecell)*r->N_tree_essential_recv_max[i]);
		}
	}
	
	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[r->mpi_num];
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_tree_essential_recv[i]==0) continue;
		MPI_Irecv(r->tree_essential_recv[i], sizeof(struct reb_treecell)* r->N_tree_essential_recv[i], MPI_CHAR, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_tree_essential_send[i]==0) continue;
		MPI_Send(r->tree_essential_send[i],  sizeof(struct reb_treecell)*r->N_tree_essential_send[i], MPI_CHAR, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_tree_essential_recv[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		for (int j=0;j<r->N_tree_essential_recv[i];j++){
			reb_tree_add_essential_node(r, &(r->tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->N_tree_essential_send[i] = 0;
		r->N_tree_essential_recv[i] = 0;
	}
}

void reb_communication_mpi_distribute_essential_tree_for_collisions(struct reb_simulation* const r){
	///////////////////////////////////////////////////////////////
	// Distribute essential tree needed for gravity and collisions
	///////////////////////////////////////////////////////////////
	
	// Distribute the number of cells to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->N_tree_essential_send, 1, MPI_INT, &(r->N_tree_essential_recv[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->N_tree_essential_recv_max[i]<r->N_tree_essential_recv[i]){
			r->N_tree_essential_recv_max[i] += 32;
			r->tree_essential_recv[i] = realloc(r->tree_essential_recv[i],sizeof(struct reb_treecell)*r->N_tree_essential_recv_max[i]);
		}
	}

	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[r->mpi_num];
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_tree_essential_recv[i]==0) continue;
		MPI_Irecv(r->tree_essential_recv[i],  sizeof(struct reb_treecell)*r->N_tree_essential_recv[i], MPI_CHAR, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_tree_essential_send[i]==0) continue;
		MPI_Send(r->tree_essential_send[i],  sizeof(struct reb_treecell)*r->N_tree_essential_send[i], MPI_CHAR, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_tree_essential_recv[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<r->mpi_num;i++){
		for (int j=0;j<r->N_tree_essential_recv[i];j++){
			reb_tree_add_essential_node(r, &(r->tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->N_tree_essential_send[i] = 0;
		r->N_tree_essential_recv[i] = 0;
	}

	//////////////////////////////////////////////////////
	// Distribute particles needed for collisiosn search 
	//////////////////////////////////////////////////////
	
	// Distribute the number of particles to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->N_particles_send, 1, MPI_INT, &(r->N_particles_recv[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->N_particles_recv_max[i]<r->N_particles_recv[i]){
			r->N_particles_recv_max[i] += 32;
			r->particles_recv[i] = realloc(r->particles_recv[i],sizeof(struct reb_particle)*r->N_particles_recv_max[i]);
		}
	}
	
	// Exchange particles via MPI.
	// Using non-blocking receive call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_particles_recv[i]==0) continue;
		MPI_Irecv(r->particles_recv[i], sizeof(struct reb_particle)* r->N_particles_recv[i], MPI_CHAR, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_particles_send[i]==0) continue;
		MPI_Send(r->particles_send[i], sizeof(struct reb_particle)* r->N_particles_send[i], MPI_CHAR, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->N_particles_recv[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// No need to add particles to tree as reference already set.
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->N_particles_send[i] = 0;
		r->N_particles_recv[i] = 0;
	}
}


#endif // MPI
