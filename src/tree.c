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

void add_particle_to_tree(int pt);

void tree_init(){
	printf("Initializing the tree\n");
	root_nx = round(boxsize_x/boxsize_min);
	root_ny = round(boxsize_y/boxsize_min);
	root_nz = round(boxsize_z/boxsize_min);
	boxsize_x = (double)root_nx * boxsize_min;
	boxsize_y = (double)root_ny * boxsize_min;
	boxsize_z = (double)root_nz * boxsize_min;
	root = malloc(root_nx*root_ny*root_nz*sizeof(struct cell*));
	for(int i=0;i<root_nx;i++){
	for(int j=0;j<root_ny;j++){
	for(int k=0;k<root_nz;k++){
		int index = (k*root_ny+j)*root_nx+i;
		root[index] = NULL;
	}
	}
	}

	for (int i=0;i<N;i++){
		add_particle_to_tree(i);
	}
}


int get_roottree_for_particle(int pt){
	struct particle p = particles[pt];
	int i = ((int)floor((p.x + boxsize_x/2.)/boxsize_min));
	int j = ((int)floor((p.y + boxsize_y/2.)/boxsize_min));
	int k = ((int)floor((p.z + boxsize_z/2.)/boxsize_min));
	int index = (k*root_ny+j)*root_nx+i;
	return index;
}

void add_particle_to_tree(int pt){
	root[get_roottree_for_particle(pt)] = add_particle(root[get_roottree_for_particle(pt)],pt,NULL,0);
}

struct cell *add_particle(struct cell *node, int n, struct cell *parent, int o){
	if (node == NULL) {
		node = calloc(1, sizeof(struct cell));
		struct particle p = particles[n];
		node->m		= p.m;
		node->mx	= p.x;
		node->my	= p.y;
		node->mz	= p.z;
		if (parent == NULL){
			node->w = boxsize_min;
			int i = ((int)floor((p.x + boxsize_x/2.)/boxsize_min));
			int j = ((int)floor((p.y + boxsize_y/2.)/boxsize_min));
			int k = ((int)floor((p.z + boxsize_z/2.)/boxsize_min));
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
		node->pt = n;
		for (int i=0; i<8; i++) node->oct[i] = NULL;
	} else if (node->pt >= 0) {
//		node.m	+= particles[n].m;
//		node.mx	=  (node.mx*node.m + particles[n].m*particles[n].x) / node.m;
//		node.my	=  (node.my*node.m + particles[n].m*particles[n].y) / node.m;
//		node.mz	=  (node.mz*node.m + particles[n].m*particles[n].z) / node.m;
		int o;
		o = get_octant(node->pt, node);
		node->oct[o] = add_particle(node->oct[o], node->pt, node, o); 
		o = get_octant(n, node);
		node->oct[o] = add_particle(node->oct[o], n, node, o);
		node->pt = -2;
		}
	else {
		node->pt--;
		int o = get_octant(n, node);
		node->oct[o] = add_particle(node->oct[o], n, node, o);
	}
	return node;
}

int get_octant(int n, struct cell *node){
	int octant = 0;
	octant = ldexp((particles[n].x > node->x) ? 0 : 1 , 0)	+ \
			ldexp((particles[n].y > node->y) ? 0 : 1 , 1) + \
			ldexp((particles[n].z > node->z) ? 0 : 1 , 2);
	return octant;
}

int isInside(struct cell *node){
	if (fabs(particles[node->pt].x-node->x) > node->w || \
		fabs(particles[node->pt].y-node->y) > node->w || \
		fabs(particles[node->pt].z-node->z) > node->w) {
		return 0;
	}
	return 1;
}

struct cell *update_tree(struct cell *node){
	int test;
	if (node == NULL) {
		return NULL;
	// Non-leaf nodes	
	} else	if (node->pt < 0) {
		for (int o = 0; o < 8; o++) {
			node->oct[o] = update_tree(node->oct[o]);
		}
		// Calculate the total mass and center of mass, and check if the node needs 
		// derefinement after updating the tree.
		node->pt = 0;
		node->m	 = 0;
		node->mx = 0;
		node->my = 0;
		node->mz = 0;
		for (int o = 0; o < 8; o++) {
			if (node->oct[o] != NULL) {
				// Calculate the total mass and the center of mass
				node->mx = (node->mx*node->m + node->oct[o]->mx* node->oct[o]->m) / (node->m + node->oct[o]->m);
				node->my = (node->my*node->m + node->oct[o]->my* node->oct[o]->m) / (node->m + node->oct[o]->m);
				node->mz = (node->mz*node->m + node->oct[o]->mz* node->oct[o]->m) / (node->m + node->oct[o]->m);
				node->m += node->oct[o]->m;
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
		if (t > 0) {
			if (node->pt == 0) {	// The node is empty.
				free(node);
				return NULL;
			} else if (node->pt == -1) { // The node becomes a leaf.
				node->pt = node->oct[test]->pt;
				free(node->oct[test]);
				node->oct[test]=NULL;
				return node;
			}
		}
		return node;
	// Leaf nodes
	} else if (t == 0) {
		return node;
	} else if (isInside(node) == 0) {
		add_particle_to_tree(node->pt);
		free(node);
		return NULL; 
	} else {
		return node;
	}
}

void tree_update(){
	for(int i=0;i<root_nx;i++){
	for(int j=0;j<root_ny;j++){
	for(int k=0;k<root_nz;k++){
		int index = (k*root_ny+j)*root_nx+i;
		root[index] = update_tree(root[index]);
	}
	}
	}
}
