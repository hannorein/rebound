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
#include "rebound.h"
#include "tree.h"
#include "boundary.h"
#include "communication_mpi.h"

void reb_communication_mpi_init(struct reb_simulation* const r, int argc, char** argv){
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&(r->mpi_num));
	MPI_Comm_rank(MPI_COMM_WORLD,&(r->mpi_id));
	
	
	// Setup MPI description of the particle structure 
	int bnum = 0;
	int blen[5];
	MPI_Aint indices[5];
	MPI_Datatype oldtypes[5];
    {
        blen[bnum] 	= 12;
        indices[bnum] 	= 0; 
        oldtypes[bnum] 	= MPI_DOUBLE;
    }
	bnum++;
    {
        struct reb_particle p;
        blen[bnum] 	= 1; 
        indices[bnum] 	= (char*)&p.c - (char*)&p; 
        oldtypes[bnum] 	= MPI_CHAR;
    }
	bnum++;
    {
        struct reb_particle p;
        blen[bnum] 	= 1; 
        indices[bnum] 	= (char*)&p.hash - (char*)&p; 
        oldtypes[bnum] 	= MPI_INT; // hash is a uint32_t, but not all MPI headers seem to have MPI_UINT32_T
    }
	bnum++;
    {
        struct reb_particle p;
        blen[bnum] 	= 2; 
        indices[bnum] 	= (char*)&p.ap - (char*)&p; 
        oldtypes[bnum] 	= MPI_CHAR;
    }
	bnum++;
    {
        blen[bnum] 	= 1; 
        indices[bnum] 	= sizeof(struct reb_particle); 
        oldtypes[bnum] 	= MPI_UB;
    }
	bnum++;
	MPI_Type_struct(bnum, blen, indices, oldtypes, &(r->mpi_particle) );
	MPI_Type_commit(&(r->mpi_particle)); 

	// Setup MPI description of the cell structure 
	struct reb_treecell c;
	bnum = 0;
    {
        blen[bnum] 	= 8; 
#ifdef QUADRUPOLE
        blen[bnum] 	+= 6;
#endif // QUADRUPOLE
        indices[bnum] 	= 0; 
        oldtypes[bnum] 	= MPI_DOUBLE;
    }
	bnum++;
    {
        blen[bnum] 	= 8; 
        indices[bnum] 	= (char*)&c.oct - (char*)&c; 
        oldtypes[bnum] 	= MPI_CHAR;
    }
	bnum++;
    {
        blen[bnum] 	= 1; 
        indices[bnum] 	= (char*)&c.pt - (char*)&c; 
        oldtypes[bnum] 	= MPI_INT;
    }
	bnum++;
    {
        blen[bnum] 	= 1; 
        indices[bnum] 	= sizeof(struct reb_treecell); 
        oldtypes[bnum] 	= MPI_UB;
    }
	bnum++;
	MPI_Type_struct(bnum, blen, indices, oldtypes, &(r->mpi_cell) );
	MPI_Type_commit(&(r->mpi_cell)); 
	
	// Prepare send/recv buffers for particles
	r->particles_send   	= calloc(r->mpi_num,sizeof(struct reb_particle*));
	r->particles_send_N 	= calloc(r->mpi_num,sizeof(int));
	r->particles_send_Nmax 	= calloc(r->mpi_num,sizeof(int));
	r->particles_recv   	= calloc(r->mpi_num,sizeof(struct reb_particle*));
	r->particles_recv_N 	= calloc(r->mpi_num,sizeof(int));
	r->particles_recv_Nmax 	= calloc(r->mpi_num,sizeof(int));

	// Prepare send/recv buffers for essential tree
	r->tree_essential_send   	= calloc(r->mpi_num,sizeof(struct reb_treecell*));
	r->tree_essential_send_N 	= calloc(r->mpi_num,sizeof(int));
	r->tree_essential_send_Nmax = calloc(r->mpi_num,sizeof(int));
	r->tree_essential_recv   	= calloc(r->mpi_num,sizeof(struct reb_treecell*));
	r->tree_essential_recv_N 	= calloc(r->mpi_num,sizeof(int));
	r->tree_essential_recv_Nmax = calloc(r->mpi_num,sizeof(int));
}

int reb_communication_mpi_rootbox_is_local(struct reb_simulation* const r, int i){
	int root_n_per_node = r->root_n/r->mpi_num;
	int proc_id = i/root_n_per_node;
	if (proc_id != r->mpi_id){
		return 0;
	}else{
		return 1;
	}
}


void reb_communication_mpi_distribute_particles(struct reb_simulation* const r){
	// Distribute the number of particles to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->particles_send_N, 1, MPI_INT, &(r->particles_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->particles_recv_Nmax[i]<r->particles_recv_N[i]){
			r->particles_recv_Nmax[i] += 32;
			r->particles_recv[i] = realloc(r->particles_recv[i],sizeof(struct reb_particle)*r->particles_recv_Nmax[i]);
		}
	}

	// Exchange particles via MPI.
	// Using non-blocking receive call.
	MPI_Request request[r->mpi_num];
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->particles_recv_N[i]==0) continue;
		MPI_Irecv(r->particles_recv[i], r->particles_recv_N[i], r->mpi_particle, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->particles_send_N[i]==0) continue;
		MPI_Send(r->particles_send[i], r->particles_send_N[i], r->mpi_particle, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->particles_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add particles to local tree
	for (int i=0;i<r->mpi_num;i++){
		for (int j=0;j<r->particles_recv_N[i];j++){
			reb_add(r,r->particles_recv[i][j]);
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->particles_send_N[i] = 0;
		r->particles_recv_N[i] = 0;
	}
}

void reb_communication_mpi_add_particle_to_send_queue(struct reb_simulation* const r, struct reb_particle pt, int proc_id){
	int send_N = r->particles_send_N[proc_id];
	while (r->particles_send_Nmax[proc_id] <= send_N){
		r->particles_send_Nmax[proc_id] += 128;
		r->particles_send[proc_id] = realloc(r->particles_send[proc_id],sizeof(struct reb_particle)*r->particles_send_Nmax[proc_id]);
	}
	r->particles_send[proc_id][send_N] = pt;
	r->particles_send_N[proc_id]++;
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
	int i = index%r->root_nx;
	int j = ((index-i)/r->root_nx)%r->root_ny;
	int k = ((index-i)/r->root_nx-j)/r->root_ny;
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
	int root_n_per_node = r->root_n/r->mpi_num;
	int root_start = proc_id*root_n_per_node;
	int root_stop  = (proc_id+1)*root_n_per_node;
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
	int nghostxcol = (r->nghostx>0?1:0);
	int nghostycol = (r->nghosty>0?1:0);
	int nghostzcol = (r->nghostz>0?1:0);
	double distance2 = r->root_size*(double)r->root_n; // A conservative estimate for the minimum distance.
	distance2 *= distance2;
	for (int gbx=-nghostxcol; gbx<=nghostxcol; gbx++){
	for (int gby=-nghostycol; gby<=nghostycol; gby++){
	for (int gbz=-nghostzcol; gbz<=nghostzcol; gbz++){
		struct reb_ghostbox gb = reb_boundary_get_ghostbox(r, gbx,gby,gbz);
		struct reb_aabb boundingbox = reb_communication_boundingbox_for_proc(r, proc_id);
		boundingbox.xmin+=gb.shiftx;
		boundingbox.xmax+=gb.shiftx;
		boundingbox.ymin+=gb.shifty;
		boundingbox.ymax+=gb.shifty;
		boundingbox.zmin+=gb.shiftz;
		boundingbox.zmax+=gb.shiftz;
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
	if (r->tree_essential_send_N[proc]>=r->tree_essential_send_Nmax[proc]){
		r->tree_essential_send_Nmax[proc] += 32;
		r->tree_essential_send[proc] = realloc(r->tree_essential_send[proc],sizeof(struct reb_treecell)*r->tree_essential_send_Nmax[proc]);
	}
	// Copy node to send buffer
	r->tree_essential_send[proc][r->tree_essential_send_N[proc]] = (*node);
	r->tree_essential_send_N[proc]++;
	if (node->pt>=0){ // Is leaf
		// Also transmit particle (Here could be another check if the particle actually overlaps with the other box)
		if (r->particles_send_N[proc]>=r->particles_send_Nmax[proc]){
			r->particles_send_Nmax[proc] += 32;
			r->particles_send[proc] = realloc(r->particles_send[proc],sizeof(struct reb_treecell)*r->particles_send_Nmax[proc]);
		}
		// Copy particle to send buffer
		r->particles_send[proc][r->particles_send_N[proc]] = r->particles[node->pt];
		// Update reference from cell to particle
		r->tree_essential_send[proc][r->tree_essential_send_N[proc]-1].pt = r->particles_send_N[proc];
		r->particles_send_N[proc]++;
	}else{		// Not a leaf. Check if we need to transfer daughters.
		double distance2 = reb_communication_distance2_of_proc_to_node(r, proc,node);
		double rp  = 2.*r->max_radius[0] + 0.86602540378443*node->w;
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
	if (r->tree_essential_send_N[proc]>=r->tree_essential_send_Nmax[proc]){
		r->tree_essential_send_Nmax[proc] += 32;
		r->tree_essential_send[proc] = realloc(r->tree_essential_send[proc],sizeof(struct reb_treecell)*r->tree_essential_send_Nmax[proc]);
	}
	// Copy node to send buffer
	r->tree_essential_send[proc][r->tree_essential_send_N[proc]] = (*node);
	r->tree_essential_send_N[proc]++;
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
		MPI_Scatter(r->tree_essential_send_N, 1, MPI_INT, &(r->tree_essential_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->tree_essential_recv_Nmax[i]<r->tree_essential_recv_N[i]){
			r->tree_essential_recv_Nmax[i] += 32;
			r->tree_essential_recv[i] = realloc(r->tree_essential_recv[i],sizeof(struct reb_treecell)*r->tree_essential_recv_Nmax[i]);
		}
	}
	
	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[r->mpi_num];
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->tree_essential_recv_N[i]==0) continue;
		MPI_Irecv(r->tree_essential_recv[i], r->tree_essential_recv_N[i], r->mpi_cell, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->tree_essential_send_N[i]==0) continue;
		MPI_Send(r->tree_essential_send[i], r->tree_essential_send_N[i], r->mpi_cell, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->tree_essential_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		for (int j=0;j<r->tree_essential_recv_N[i];j++){
			reb_tree_add_essential_node(r, &(r->tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->tree_essential_send_N[i] = 0;
		r->tree_essential_recv_N[i] = 0;
	}
}

void reb_communication_mpi_distribute_essential_tree_for_collisions(struct reb_simulation* const r){
	///////////////////////////////////////////////////////////////
	// Distribute essential tree needed for gravity and collisions
	///////////////////////////////////////////////////////////////
	
	// Distribute the number of cells to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->tree_essential_send_N, 1, MPI_INT, &(r->tree_essential_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming tree_essential
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->tree_essential_recv_Nmax[i]<r->tree_essential_recv_N[i]){
			r->tree_essential_recv_Nmax[i] += 32;
			r->tree_essential_recv[i] = realloc(r->tree_essential_recv[i],sizeof(struct reb_treecell)*r->tree_essential_recv_Nmax[i]);
		}
	}

	// Exchange tree_essential via MPI.
	// Using non-blocking receive call.
	MPI_Request request[r->mpi_num];
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->tree_essential_recv_N[i]==0) continue;
		MPI_Irecv(r->tree_essential_recv[i], r->tree_essential_recv_N[i], r->mpi_cell, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->tree_essential_send_N[i]==0) continue;
		MPI_Send(r->tree_essential_send[i], r->tree_essential_send_N[i], r->mpi_cell, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all tree_essential to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->tree_essential_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// Add tree_essential to local tree
	for (int i=0;i<r->mpi_num;i++){
		for (int j=0;j<r->tree_essential_recv_N[i];j++){
			reb_tree_add_essential_node(r, &(r->tree_essential_recv[i][j]));
		}
	}
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->tree_essential_send_N[i] = 0;
		r->tree_essential_recv_N[i] = 0;
	}

	//////////////////////////////////////////////////////
	// Distribute particles needed for collisiosn search 
	//////////////////////////////////////////////////////
	
	// Distribute the number of particles to be transferred.
	for (int i=0;i<r->mpi_num;i++){
		MPI_Scatter(r->particles_send_N, 1, MPI_INT, &(r->particles_recv_N[i]), 1, MPI_INT, i, MPI_COMM_WORLD);
	}
	// Allocate memory for incoming particles
	for (int i=0;i<r->mpi_num;i++){
		if  (i==r->mpi_id) continue;
		while (r->particles_recv_Nmax[i]<r->particles_recv_N[i]){
			r->particles_recv_Nmax[i] += 32;
			r->particles_recv[i] = realloc(r->particles_recv[i],sizeof(struct reb_particle)*r->particles_recv_Nmax[i]);
		}
	}
	
	// Exchange particles via MPI.
	// Using non-blocking receive call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->particles_recv_N[i]==0) continue;
		MPI_Irecv(r->particles_recv[i], r->particles_recv_N[i], r->mpi_particle, i, i*r->mpi_num+r->mpi_id, MPI_COMM_WORLD, &(request[i]));
	}
	// Using blocking send call.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->particles_send_N[i]==0) continue;
		MPI_Send(r->particles_send[i], r->particles_send_N[i], r->mpi_particle, i, r->mpi_id*r->mpi_num+i, MPI_COMM_WORLD);
	}
	// Wait for all particles to be received.
	for (int i=0;i<r->mpi_num;i++){
		if (i==r->mpi_id) continue;
		if (r->particles_recv_N[i]==0) continue;
		MPI_Status status;
		MPI_Wait(&(request[i]), &status);
	}
	// No need to add particles to tree as reference already set.
	// Bring everybody into sync, clean up. 
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0;i<r->mpi_num;i++){
		r->particles_send_N[i] = 0;
		r->particles_recv_N[i] = 0;
	}
}


#endif // MPI
