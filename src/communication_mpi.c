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
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"
#include "collisions.h"
#include "communication_mpi.h"

#include "mpi.h"
MPI_Datatype mpi_particle;
#ifdef TREE
MPI_Datatype mpi_cell;
#endif
MPI_Status stat; 
int mpi_num;
int mpi_id;

struct particle** 	particles_send;
int* 			particles_send_N;
int* 			particles_send_Nmax;
struct particle** 	particles_recv;
int* 			particles_recv_N;
int* 			particles_recv_Nmax;

#ifdef TREE
struct cell** 	tree_essential_send;
int* 		tree_essential_send_N;
int* 		tree_essential_send_Nmax;
struct cell** 	tree_essential_recv;
int* 		tree_essential_recv_N;
int* 		tree_essential_recv_Nmax;
#endif


void communication_mpi_init(int argc, char** argv){
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_num);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
	
	
	// Setup MPI description of the particle structure 
	int bnum = 0;
	int blen[4];
	MPI_Aint indices[4];
	MPI_Datatype oldtypes[4];
	blen[bnum] 	= 10;
#ifndef COLLISIONS_NONE
	blen[bnum] 	+= 2; 
#endif //COLLISIONS_NONE
	indices[bnum] 	= 0; 
	oldtypes[bnum] 	= MPI_DOUBLE;
	bnum++;
#ifdef TREE
	struct particle p;
	blen[bnum] 	= 1; 
	indices[bnum] 	= (char*)&p.c - (char*)&p; 
	oldtypes[bnum] 	= MPI_CHAR;
	bnum++;
#endif // TREE
	blen[bnum] 	= 1; 
	indices[bnum] 	= sizeof(struct particle); 
	oldtypes[bnum] 	= MPI_UB;
	bnum++;
	MPI_Type_struct(bnum, blen, indices, oldtypes, &mpi_particle );
	MPI_Type_commit(&mpi_particle); 

#ifdef TREE
	// Setup MPI description of the cell structure 
	struct cell c;
	bnum = 0;
	blen[bnum] 	= 4; 
#ifdef GRAVITY_TREE
	blen[bnum] 	+= 4;
#ifdef QUADRUPOLE
	blen[bnum] 	+= 6;
#endif // QUADRUPOLE
#endif //GRAVITY_TREE
	indices[bnum] 	= 0; 
	oldtypes[bnum] 	= MPI_DOUBLE;
	bnum++;
	blen[bnum] 	= 8; 
	indices[bnum] 	= (char*)&c.oct - (char*)&c; 
	oldtypes[bnum] 	= MPI_CHAR;
	bnum++;
	blen[bnum] 	= 1; 
	indices[bnum] 	= (char*)&c.pt - (char*)&c; 
	oldtypes[bnum] 	= MPI_INT;
	bnum++;
	blen[bnum] 	= 1; 
	indices[bnum] 	= sizeof(struct cell); 
	oldtypes[bnum] 	= MPI_UB;
	bnum++;
	MPI_Type_struct(bnum, blen, indices, oldtypes, &mpi_cell );
	MPI_Type_commit(&mpi_cell); 
#endif // TREE 
	
	// Prepare send/recv buffers for particles
	particles_send   	= calloc(mpi_num,sizeof(struct particle*));
	particles_send_N 	= calloc(mpi_num,sizeof(int));
	particles_send_Nmax 	= calloc(mpi_num,sizeof(int));
	particles_recv   	= calloc(mpi_num,sizeof(struct particle*));
	particles_recv_N 	= calloc(mpi_num,sizeof(int));
	particles_recv_Nmax 	= calloc(mpi_num,sizeof(int));

#ifdef TREE
	// Prepare send/recv buffers for essential tree
	tree_essential_send   	= calloc(mpi_num,sizeof(struct cell*));
	tree_essential_send_N 	= calloc(mpi_num,sizeof(int));
	tree_essential_send_Nmax= calloc(mpi_num,sizeof(int));
	tree_essential_recv   	= calloc(mpi_num,sizeof(struct cell*));
	tree_essential_recv_N 	= calloc(mpi_num,sizeof(int));
	tree_essential_recv_Nmax= calloc(mpi_num,sizeof(int));
#endif
}

int communication_mpi_rootbox_is_local(int i){
	int root_n_per_node = root_n/mpi_num;
	int proc_id = i/root_n_per_node;
	if (proc_id != mpi_id){
		return 0;
	}else{
		return 1;
	}
}


void communication_mpi_distribute_particles(){
	// Distribute the number of particles to be transferred.
	for (int i=0;i<mpi_num;i++){
		MPI_Scatter(particles_send_N, 1, MPI_INT, &(particles_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<mpi_num;i++){
		if  (i==mpi_id) continue;
		while (particles_recv_Nmax[i]<particles_recv_N[i]){
			particles_recv_Nmax[i] += 32;
			particles_recv[i] = realloc(particles_recv[i],sizeof(struct particle)*particles_recv_Nmax[i]);
		}
	}

	// Exchange particles via MPI.
	// Using non-blocking receive call.
	MPI_Request request[mpi_num];
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
		MPI_Irecv(particles_recv[i], particles_recv_N[i], mpi_particle, i, i*mpi_num+mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_send_N[i]==0) continue;
		MPI_Send(particles_send[i], particles_send_N[i], mpi_particle, i, mpi_id*mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add particles to local tree
	for (int i=0;i<mpi_num;i++){
		for (int j=0;j<particles_recv_N[i];j++){
			particles_add(particles_recv[i][j]);
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		particles_send_N[i] = 0;
		particles_recv_N[i] = 0;
	}
}

void communication_mpi_add_particle_to_send_queue(struct particle pt, int proc_id){
	int send_N = particles_send_N[proc_id];
	while (particles_send_Nmax[proc_id] <= send_N){
		particles_send_Nmax[proc_id] += 128;
		particles_send[proc_id] = realloc(particles_send[proc_id],sizeof(struct particle)*particles_send_Nmax[proc_id]);
	}
	particles_send[proc_id][send_N] = pt;
	particles_send_N[proc_id]++;
}

#ifdef TREE

/** 
 * This is the data structure for an axis aligned bounding box.
 */
struct aabb{
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
};

struct aabb communication_boundingbox_for_root(int index){
	int i = index%root_nx;
	int j = ((index-i)/root_nx)%root_ny;
	int k = ((index-i)/root_nx-j)/root_ny;
	struct aabb boundingbox;
	boundingbox.xmin = -boxsize_x/2.+boxsize*(double)i;
	boundingbox.ymin = -boxsize_y/2.+boxsize*(double)j;
	boundingbox.zmin = -boxsize_z/2.+boxsize*(double)k;
	boundingbox.xmax = -boxsize_x/2.+boxsize*(double)(i+1);
	boundingbox.ymax = -boxsize_y/2.+boxsize*(double)(j+1);
	boundingbox.zmax = -boxsize_z/2.+boxsize*(double)(k+1);
	return boundingbox;
}

struct aabb communication_boundingbox_for_proc(int proc_id){
	int root_n_per_node = root_n/mpi_num;
	int root_start = proc_id*root_n_per_node;
	int root_stop  = (proc_id+1)*root_n_per_node;
	struct aabb boundingbox = communication_boundingbox_for_root(root_start);
	for (int i=root_start+1;i<root_stop;i++){
		struct aabb boundingbox2 = communication_boundingbox_for_root(i);
		if (boundingbox.xmin > boundingbox2.xmin) boundingbox.xmin = boundingbox2.xmin;
		if (boundingbox.ymin > boundingbox2.ymin) boundingbox.ymin = boundingbox2.ymin;
		if (boundingbox.zmin > boundingbox2.zmin) boundingbox.zmin = boundingbox2.zmin;
		if (boundingbox.xmax < boundingbox2.xmax) boundingbox.xmax = boundingbox2.xmax;
		if (boundingbox.ymax < boundingbox2.ymax) boundingbox.ymax = boundingbox2.ymax;
		if (boundingbox.zmax < boundingbox2.zmax) boundingbox.zmax = boundingbox2.zmax;
	}
	return boundingbox;
}

double communication_distance2_of_aabb_to_cell(struct aabb bb, struct cell* node){
	double distancex = fabs(node->x - (bb.xmin+bb.xmax)/2.)  -  (node->w + bb.xmax-bb.xmin)/2.;
	double distancey = fabs(node->y - (bb.ymin+bb.ymax)/2.)  -  (node->w + bb.ymax-bb.ymin)/2.;
	double distancez = fabs(node->z - (bb.zmin+bb.zmax)/2.)  -  (node->w + bb.zmax-bb.zmin)/2.;
	if (distancex<0) distancex =0;
	if (distancey<0) distancey =0;
	if (distancez<0) distancez =0;
	return distancex*distancex + distancey*distancey + distancez*distancez;
}

double communication_distance2_of_proc_to_node(int proc_id, struct cell* node){
	int nghostxcol = (nghostx>0?1:0);
	int nghostycol = (nghosty>0?1:0);
	int nghostzcol = (nghostz>0?1:0);
	double distance2 = boxsize*(double)root_n; // A conservative estimate for the minimum distance.
	distance2 *= distance2;
	for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
	for (int gby=-nghostycol; gby<=nghostycol; gby++){
	for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
		struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,gbz);
		struct aabb boundingbox = communication_boundingbox_for_proc(proc_id);
		boundingbox.xmin+=gb.shiftx;
		boundingbox.xmax+=gb.shiftx;
		boundingbox.ymin+=gb.shifty;
		boundingbox.ymax+=gb.shifty;
		boundingbox.zmin+=gb.shiftz;
		boundingbox.zmax+=gb.shiftz;
		// calculate distance
		double distance2new = communication_distance2_of_aabb_to_cell(boundingbox,node);
		if (distance2 > distance2new) distance2 = distance2new;
	}
	}
	}
	return distance2;
}

void communication_mpi_prepare_essential_cell_for_collisions_for_proc(struct cell* node, int proc){
	// Add essential cell to tree_essential_send
	if (tree_essential_send_N[proc]>=tree_essential_send_Nmax[proc]){
		tree_essential_send_Nmax[proc] += 32;
		tree_essential_send[proc] = realloc(tree_essential_send[proc],sizeof(struct cell)*tree_essential_send_Nmax[proc]);
	}
	// Copy node to send buffer
	tree_essential_send[proc][tree_essential_send_N[proc]] = (*node);
	tree_essential_send_N[proc]++;
	if (node->pt>=0){ // Is leaf
		// Also transmit particle (Here could be another check if the particle actually overlaps with the other box)
		if (particles_send_N[proc]>=particles_send_Nmax[proc]){
			particles_send_Nmax[proc] += 32;
			particles_send[proc] = realloc(particles_send[proc],sizeof(struct cell)*particles_send_Nmax[proc]);
		}
		// Copy particle to send buffer
		particles_send[proc][particles_send_N[proc]] = particles[node->pt];
		// Update reference from cell to particle
		tree_essential_send[proc][tree_essential_send_N[proc]-1].pt = particles_send_N[proc];
		particles_send_N[proc]++;
	}else{		// Not a leaf. Check if we need to transfer daughters.
		double distance2 = communication_distance2_of_proc_to_node(proc,node);
		double rp  = 2.*collisions_max_r + 0.86602540378443*node->w;
		if (distance2 < rp*rp ){
			for (int o=0;o<8;o++){
				struct cell* d = node->oct[o];
				if (d==NULL) continue;
				communication_mpi_prepare_essential_cell_for_collisions_for_proc(d,proc);
			}
		}
	}
}
void communication_mpi_prepare_essential_tree_for_collisions(struct cell* root){
	if (root==NULL) return;
	// Find out which cells are needed by every other node
	for (int i=0; i<mpi_num; i++){
		if (i==mpi_id) continue;
		communication_mpi_prepare_essential_cell_for_collisions_for_proc(root,i);	
	}
}


#ifdef GRAVITY_TREE
extern double opening_angle2;
void communication_mpi_prepare_essential_cell_for_gravity_for_proc(struct cell* node, int proc){
	// Add essential cell to tree_essential_send
	if (tree_essential_send_N[proc]>=tree_essential_send_Nmax[proc]){
		tree_essential_send_Nmax[proc] += 32;
		tree_essential_send[proc] = realloc(tree_essential_send[proc],sizeof(struct cell)*tree_essential_send_Nmax[proc]);
	}
	// Copy node to send buffer
	tree_essential_send[proc][tree_essential_send_N[proc]] = (*node);
	tree_essential_send_N[proc]++;
	if (node->pt<0){		// Not a leaf. Check if we need to transfer daughters.
		double width = node->w;
		double distance2 = communication_distance2_of_proc_to_node(proc,node);
		if ( width*width > opening_angle2*distance2) {
			for (int o=0;o<8;o++){
				struct cell* d = node->oct[o];
				if (d!=NULL){
					communication_mpi_prepare_essential_cell_for_gravity_for_proc(d,proc);
				}
			}
		}
	}
}

void communication_mpi_prepare_essential_tree_for_gravity(struct cell* root){
	if (root==NULL) return;
	// Find out which cells are needed by every other node
	for (int i=0; i<mpi_num; i++){
		if (i==mpi_id) continue;
		communication_mpi_prepare_essential_cell_for_gravity_for_proc(root,i);	
	}
}

void communication_mpi_distribute_essential_tree_for_gravity(){
	///////////////////////////////////////////////////////////////
	// Distribute essential tree needed for gravity and collisions
	///////////////////////////////////////////////////////////////
	
	// Distribute the number of cells to be transferred.
	for (int i=0;i<mpi_num;i++){
		MPI_Scatter(tree_essential_send_N, 1, MPI_INT, &(tree_essential_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<mpi_num;i++){
		if  (i==mpi_id) continue;
		while (tree_essential_recv_Nmax[i]<tree_essential_recv_N[i]){
			tree_essential_recv_Nmax[i] += 32;
			tree_essential_recv[i] = realloc(tree_essential_recv[i],sizeof(struct cell)*tree_essential_recv_Nmax[i]);
		}
	}
	
	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[mpi_num];
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_recv_N[i]==0) continue;
		MPI_Irecv(tree_essential_recv[i], tree_essential_recv_N[i], mpi_cell, i, i*mpi_num+mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_send_N[i]==0) continue;
		MPI_Send(tree_essential_send[i], tree_essential_send_N[i], mpi_cell, i, mpi_id*mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		for (int j=0;j<tree_essential_recv_N[i];j++){
			tree_add_essential_node(&(tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		tree_essential_send_N[i] = 0;
		tree_essential_recv_N[i] = 0;
	}
}
#endif // GRAVITY_TREE

void communication_mpi_distribute_essential_tree_for_collisions(){
	///////////////////////////////////////////////////////////////
	// Distribute essential tree needed for gravity and collisions
	///////////////////////////////////////////////////////////////
	
	// Distribute the number of cells to be transferred.
	for (int i=0;i<mpi_num;i++){
		MPI_Scatter(tree_essential_send_N, 1, MPI_INT, &(tree_essential_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<mpi_num;i++){
		if  (i==mpi_id) continue;
		while (tree_essential_recv_Nmax[i]<tree_essential_recv_N[i]){
			tree_essential_recv_Nmax[i] += 32;
			tree_essential_recv[i] = realloc(tree_essential_recv[i],sizeof(struct cell)*tree_essential_recv_Nmax[i]);
		}
	}

	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[mpi_num];
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_recv_N[i]==0) continue;
		MPI_Irecv(tree_essential_recv[i], tree_essential_recv_N[i], mpi_cell, i, i*mpi_num+mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_send_N[i]==0) continue;
		MPI_Send(tree_essential_send[i], tree_essential_send_N[i], mpi_cell, i, mpi_id*mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (tree_essential_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<mpi_num;i++){
		for (int j=0;j<tree_essential_recv_N[i];j++){
			tree_add_essential_node(&(tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		tree_essential_send_N[i] = 0;
		tree_essential_recv_N[i] = 0;
	}

	//////////////////////////////////////////////////////
	// Distribute particles needed for collisiosn search 
	//////////////////////////////////////////////////////
	
	// Distribute the number of particles to be transferred.
	for (int i=0;i<mpi_num;i++){
		MPI_Scatter(particles_send_N, 1, MPI_INT, &(particles_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<mpi_num;i++){
		if  (i==mpi_id) continue;
		while (particles_recv_Nmax[i]<particles_recv_N[i]){
			particles_recv_Nmax[i] += 32;
			particles_recv[i] = realloc(particles_recv[i],sizeof(struct particle)*particles_recv_Nmax[i]);
		}
	}
	
	// Exchange particles via MPI.
	// Using non-blocking receive call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
		MPI_Irecv(particles_recv[i], particles_recv_N[i], mpi_particle, i, i*mpi_num+mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_send_N[i]==0) continue;
		MPI_Send(particles_send[i], particles_send_N[i], mpi_particle, i, mpi_id*mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<mpi_num;i++){
		if (i==mpi_id) continue;
		if (particles_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// No need to add particles to tree as reference already set.
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<mpi_num;i++){
		particles_send_N[i] = 0;
		particles_recv_N[i] = 0;
	}
}
#endif


#endif // MPI
