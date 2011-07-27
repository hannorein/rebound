#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"
#include "communication_mpi.h"

#ifdef TREE

struct cell** tree_root;

int tree_get_octant_for_particle_in_cell(int pt, struct cell *node);
struct cell *tree_add_particle_to_cell(struct cell *node, int pt, struct cell *parent, int o);

void tree_add_particle_to_tree(int pt){
	if (tree_root==NULL){
		tree_root = calloc(root_nx*root_ny*root_nz,sizeof(struct cell*));
	}
	struct particle p = particles[pt];
	int root_index = particles_get_rootbox_for_particle(p);
	tree_root[root_index] = tree_add_particle_to_cell(tree_root[root_index],pt,NULL,0);
}

struct cell *tree_add_particle_to_cell(struct cell *node, int pt, struct cell *parent, int o){
	if (node == NULL) {
		node = calloc(1, sizeof(struct cell));
		struct particle p = particles[pt];
		if (parent == NULL){
			node->w = boxsize;
			int i = ((int)floor((p.x + boxsize_x/2.)/boxsize))%root_nx;
			int j = ((int)floor((p.y + boxsize_y/2.)/boxsize))%root_ny;
			int k = ((int)floor((p.z + boxsize_z/2.)/boxsize))%root_nz;
			node->x = -boxsize_x/2.+boxsize*(0.5+(double)i);
			node->y = -boxsize_y/2.+boxsize*(0.5+(double)j);
			node->z = -boxsize_z/2.+boxsize*(0.5+(double)k);
		}else{
			node->w 	= parent->w/2.;
			node->x 	= parent->x + node->w/2.*((o>>0)%2==0?1.:-1);
			node->y 	= parent->y + node->w/2.*((o>>1)%2==0?1.:-1);
			node->z 	= parent->z + node->w/2.*((o>>2)%2==0?1.:-1);
		}

		// Double usages: in a leaf node, it stores the index of a particle; in a non-
		// leaf node, it equals to (-1)*Total Number of particles within that cell.
		node->pt = pt;
		particles[pt].c = node;
		for (int i=0; i<8; i++){
			node->oct[i] = NULL;
		}
		return node;
	}

	// Modify the total mass and center of mass of a cell
	if (node->pt >= 0) {
		int o = tree_get_octant_for_particle_in_cell(node->pt, node);
		node->oct[o] = tree_add_particle_to_cell(node->oct[o], node->pt, node, o); 
		o = tree_get_octant_for_particle_in_cell(pt, node);
		node->oct[o] = tree_add_particle_to_cell(node->oct[o], pt, node, o);
		node->pt = -2;
	}else{
		node->pt--;
		int o = tree_get_octant_for_particle_in_cell(pt, node);
		node->oct[o] = tree_add_particle_to_cell(node->oct[o], pt, node, o);
	}
	return node;
}

int tree_get_octant_for_particle_in_cell(int pt, struct cell *node){
	int octant = 0;
	struct particle p = particles[pt];
	if (p.x < node->x) octant+=1;
	if (p.y < node->y) octant+=2;
	if (p.z < node->z) octant+=4;
	return octant;
}

int tree_particle_is_inside_cell(struct cell *node){
	if (fabs(particles[node->pt].x-node->x) > node->w/2. || \
		fabs(particles[node->pt].y-node->y) > node->w/2. || \
		fabs(particles[node->pt].z-node->z) > node->w/2.) {
		return 0;
	}
	return 1;
}

struct cell *tree_update_cell(struct cell *node){
	int test = -1;
	if (node == NULL) {
		return NULL;
	}
	// Non-leaf nodes	
	if (node->pt < 0) {
		for (int o=0; o<8; o++) {
			node->oct[o] = tree_update_cell(node->oct[o]);
		}
		node->pt = 0;
		for (int o=0; o<8; o++) {
			struct cell *d = node->oct[o];
			if (d != NULL) {
				// Update node->pt
				if (d->pt >= 0) {	// The child is a leaf
					node->pt--;
					test = o;
				}else{				// The child cell contains several particles
					node->pt += d->pt;
				}
			}		
		}
		// Check if the node require derefinement.
		if (node->pt == 0) {	// The node is empty.
			free(node);
			return NULL;
		} else if (node->pt == -1) { // The node becomes a leaf.
			node->pt = node->oct[test]->pt;
			particles[node->pt].c = node;
			free(node->oct[test]);
			node->oct[test]=NULL;
			return node;
		}
		return node;
	} 
	// Leaf nodes
	if (tree_particle_is_inside_cell(node) == 0) {
		int oldpos = node->pt;
		struct particle reinsertme = particles[oldpos];
		N--;
		particles[oldpos] = particles[N];
		particles[oldpos].c->pt = oldpos;
		particles_add(reinsertme);
		free(node);
		return NULL; 
	} else {
		particles[node->pt].c = node;
		return node;
	}
}

#ifdef GRAVITY_TREE
void tree_update_gravity_data_in_cell(struct cell *node){
#ifdef QUADRUPOLE
	node->mxx = 0;
	node->mxy = 0;
	node->mxz = 0;
	node->myy = 0;
	node->myz = 0;
	node->mzz = 0;
#endif // QUADRUPOLE
	if (node->pt < 0) {
		// Non-leaf nodes	
		node->m  = 0;
		node->mx = 0;
		node->my = 0;
		node->mz = 0;
		for (int o=0; o<8; o++) {
			struct cell* d = node->oct[o];
			if (d!=NULL){
				tree_update_gravity_data_in_cell(d);
				// Calculate the total mass and the center of mass
				double d_m = d->m;
				node->mx += d->mx*d_m;
				node->my += d->my*d_m;
				node->mz += d->mz*d_m;
				node->m  += d_m;
			}
		}
		double m_tot = node->m;
		if (m_tot>0){
			node->mx /= m_tot;
			node->my /= m_tot;
			node->mz /= m_tot;
		}
#ifdef QUADRUPOLE
		for (int o=0; o<8; o++) {
			struct cell* d = node->oct[o];
			if (d!=NULL){
				// Ref: Hernquist, L., 1987, APJS
				double d_m = d->m;
				double qx  = d->mx - node->mx;
				double qy  = d->my - node->my;
				double qz  = d->mz - node->mz;
				double qr2 = qx*qx + qy*qy + qz*qz;
				node->mxx += d->mxx + d_m*(3.*qx*qx - qr2);
				node->mxy += d->mxy + d_m*3.*qx*qy;
				node->mxz += d->mxz + d_m*3.*qx*qz;
				node->myy += d->myy + d_m*(3.*qy*qy - qr2);
				node->myz += d->myz + d_m*3.*qy*qz;
			}
		}
		node->mzz = -node->mxx -node->myy;
#endif // QUADRUPOLE
	}else{ 
		// Leaf nodes
		struct particle p = particles[node->pt];
		node->m = p.m;
		node->mx = p.x;
		node->my = p.y;
		node->mz = p.z;
	}
}

void tree_update_gravity_data(){
	for(int i=0;i<root_n;i++){
#ifdef MPI
		if (communication_mpi_rootbox_is_local(i)==1){
#endif // MPI
			if (tree_root[i]!=NULL){
				tree_update_gravity_data_in_cell(tree_root[i]);
			}
#ifdef MPI
		}
#endif // MPI
	}
}
#endif // GRAVITY_TREE

void tree_update(){
	if (tree_root==NULL){
		tree_root = calloc(root_nx*root_ny*root_nz,sizeof(struct cell*));
	}
	for(int i=0;i<root_n;i++){
#ifdef MPI
		if (communication_mpi_rootbox_is_local(i)==1){
#endif // MPI
			tree_root[i] = tree_update_cell(tree_root[i]);
#ifdef MPI
		}
#endif // MPI
	}
}



#ifdef MPI
int particles_get_rootbox_for_node(struct cell* node){
	int i = ((int)floor((node->x + boxsize_x/2.)/boxsize)+root_nx)%root_nx;
	int j = ((int)floor((node->y + boxsize_y/2.)/boxsize)+root_ny)%root_ny;
	int k = ((int)floor((node->z + boxsize_z/2.)/boxsize)+root_nz)%root_nz;
	int index = (k*root_ny+j)*root_nx+i;
	return index;
}

int tree_get_octant_for_cell_in_cell(struct cell* nnode, struct cell *node){
	int octant = 0;
	if (nnode->x < node->x) octant+=1;
	if (nnode->y < node->y) octant+=2;
	if (nnode->z < node->z) octant+=4;
	return octant;
}

void tree_add_essential_node_to_node(struct cell* nnode, struct cell* node){
	int o = tree_get_octant_for_cell_in_cell(nnode, node);
	if (node->oct[o]==NULL){
		node->oct[o] = nnode;
	}else{
		tree_add_essential_node_to_node(nnode, node->oct[o]);
	}
}

void tree_add_essential_node(struct cell* node){
	// Add essential node to appropriate parent.
	for (int o=0;o<8;o++){
		node->oct[o] = NULL;	
	}
	int index = particles_get_rootbox_for_node(node);
	if (tree_root[index]==NULL){
		tree_root[index] = node;
	}else{
		tree_add_essential_node_to_node(node, tree_root[index]);
	}
}
void tree_prepare_essential_tree_for_gravity(){
	for(int i=0;i<root_n;i++){
		if (communication_mpi_rootbox_is_local(i)==1){
			communication_mpi_prepare_essential_tree_for_gravity(tree_root[i]);
		}else{
			// Delete essential tree reference. 
			// Tree itself is saved in tree_essential_recv[][] and
			// will be overwritten the next timestep.
			tree_root[i] = NULL;
		}
	}
}
void tree_prepare_essential_tree_for_collisions(){
	for(int i=0;i<root_n;i++){
		if (communication_mpi_rootbox_is_local(i)==1){
			communication_mpi_prepare_essential_tree_for_collisions(tree_root[i]);
		}else{
			// Delete essential tree reference. 
			// Tree itself is saved in tree_essential_recv[][] and
			// will be overwritten the next timestep.
			tree_root[i] = NULL;
		}
	}
}
#endif // MPI

#endif
