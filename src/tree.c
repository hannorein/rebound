#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"

struct cell** root;
int root_nx;
int root_ny;
int root_nz;

void tree_add_particle_to_tree(int pt);
int tree_get_root_for_particle(int pt);
int tree_get_octant_for_particle_in_cell(int pt, struct cell *node);
struct cell *tree_add_particle_to_cell(struct cell *node, int pt, struct cell *parent, int o);

void tree_init(){
	root_nx = round(boxsize_x/boxsize_min);
	root_ny = round(boxsize_y/boxsize_min);
	root_nz = round(boxsize_z/boxsize_min);
	boxsize_x = (double)root_nx * boxsize_min;
	boxsize_y = (double)root_ny * boxsize_min;
	boxsize_z = (double)root_nz * boxsize_min;
	root = calloc(root_nx*root_ny*root_nz,sizeof(struct cell*));

	for (int i=0;i<N;i++){
		tree_add_particle_to_tree(i);
	}
}


int tree_get_root_for_particle(int pt){
	struct particle p = particles[pt];
	int i = ((int)floor((p.x + boxsize_x/2.)/boxsize_min)+root_nx)%root_nx;
	int j = ((int)floor((p.y + boxsize_y/2.)/boxsize_min)+root_ny)%root_ny;
	int k = ((int)floor((p.z + boxsize_z/2.)/boxsize_min)+root_nz)%root_nz;
	int index = (k*root_ny+j)*root_nx+i;
	return index;
}

void tree_add_particle_to_tree(int pt){
	int root_index = tree_get_root_for_particle(pt);
	root[root_index] = tree_add_particle_to_cell(root[root_index],pt,NULL,0);
}

struct cell *tree_add_particle_to_cell(struct cell *node, int pt, struct cell *parent, int o){
	if (node == NULL) {
		node = calloc(1, sizeof(struct cell));
		struct particle p = particles[pt];
#ifdef GRAVITY_TREE
		node->m		= p.m;
		node->mx	= p.x;
		node->my	= p.y;
		node->mz	= p.z;
#endif // GRAVITY_TREE
		if (parent == NULL){
			node->w = boxsize_min;
			int i = ((int)floor((p.x + boxsize_x/2.)/boxsize_min))%root_nx;
			int j = ((int)floor((p.y + boxsize_y/2.)/boxsize_min))%root_ny;
			int k = ((int)floor((p.z + boxsize_z/2.)/boxsize_min))%root_nz;
			node->x = -boxsize_x/2.+boxsize_min*(0.5+(double)i);
			node->y = -boxsize_y/2.+boxsize_min*(0.5+(double)j);
			node->z = -boxsize_z/2.+boxsize_min*(0.5+(double)k);
		}else{
			node->w 	= parent->w/2.;
			node->x 	= parent->x + node->w/2.*((o>>0)%2==0?1.:-1);
			node->y 	= parent->y + node->w/2.*((o>>1)%2==0?1.:-1);
			node->z 	= parent->z + node->w/2.*((o>>2)%2==0?1.:-1);
		}

		// Double usages: in a leaf node, it stores the index of a particle; in a non-
		// leaf node, it equals to (-1)*Total Number of particles within that cell.
		node->pt = pt;
		for (int i=0; i<8; i++){
			node->oct[i] = NULL;
		}
		return node;
	}
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
	if (fabs(particles[node->pt].x-node->x) > node->w || \
		fabs(particles[node->pt].y-node->y) > node->w || \
		fabs(particles[node->pt].z-node->z) > node->w) {
		return 0;
	}
	return 1;
}

struct cell *tree_update_cell(struct cell *node){
	int test;
	if (node == NULL) {
		return NULL;
	}
	// Non-leaf nodes	
	if (node->pt < 0) {
		for (int o = 0; o < 8; o++) {
			node->oct[o] = tree_update_cell(node->oct[o]);
		}
		// Calculate the total mass and center of mass, and check if the node needs 
		// derefinement after updating the tree.
		node->pt = 0;
#ifdef GRAVITY_TREE
		node->m	 = 0;
		node->mx = 0;
		node->my = 0;
		node->mz = 0;
#endif // GRAVITY_TREE
		for (int o = 0; o < 8; o++) {
			if (node->oct[o] != NULL) {
				// Calculate the total mass and the center of mass
#ifdef GRAVITY_TREE
				node->mx = (node->mx*node->m + node->oct[o]->mx* node->oct[o]->m) / (node->m + node->oct[o]->m);
				node->my = (node->my*node->m + node->oct[o]->my* node->oct[o]->m) / (node->m + node->oct[o]->m);
				node->mz = (node->mz*node->m + node->oct[o]->mz* node->oct[o]->m) / (node->m + node->oct[o]->m);
				node->m += node->oct[o]->m;
#endif // GRAVITY_TREE
				// Update node->pt
				if (node->oct[o]->pt >= 0) {	// The child is a leaf
					node->pt--;
					test = o;
				}
				else {	// The child cell contains several particles
					node->pt += node->oct[o]->pt;
				}
			}		
		}
		// Check if the node require derefinement.
		if (node->pt == 0) {	// The node is empty.
			free(node);
			return NULL;
		} else if (node->pt == -1) { // The node becomes a leaf.
			node->pt = node->oct[test]->pt;
			free(node->oct[test]);
			node->oct[test]=NULL;
			return node;
		}
		return node;
	} 
	// Leaf nodes
	if (tree_particle_is_inside_cell(node) == 0) {
		tree_add_particle_to_tree(node->pt);
		free(node);
		return NULL; 
	} else {
#ifdef GRAVITY_TREE
		struct particle p = particles[node->pt];
		node->mx	= p.x;
		node->my	= p.y;
		node->mz	= p.z;
#endif // GRAVITY_TREE
		return node;
	}
}

void tree_update(){
	if (root==NULL){
		tree_init();
	}
	check_boundaries();
	for(int i=0;i<root_nx;i++){
	for(int j=0;j<root_ny;j++){
	for(int k=0;k<root_nz;k++){
		int index = (k*root_ny+j)*root_nx+i;
		root[index] = tree_update_cell(root[index]);
	}
	}
	}
}
