#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"

struct cell* root;

void init_tree(){
	printf("Initializing the tree\n");
	root = NULL;
	for (int i=0;i<N;i++){
		root = add_particle(root, i, NULL, -1);
	}
}

struct cell *add_particle(struct cell *node, int n, struct cell *parent, int o){
	if (node == NULL) {
		node = calloc(1, sizeof(struct cell));
		node->m	= particles[n].m;
		node->mx[0]	= particles[n].x;
		node->mx[1]	= particles[n].y;
		node->mx[2]	= particles[n].z;
		if (parent == NULL){
			for (int i=0; i<3; i++) node->x[i] = 0;
			node->w = boxsize;
		} else {
			node->w = parent->w/2.;
			for (int i=0; i<3; i++) {
				node->x[i] = parent->x[i] + node->w/2.*((o>>i)%2==0?1.:-1);
			}
		/*	node->x[0] = parent->x[0] + node->w/2.*(o%2==0?1.:-1);
			node->x[1] = parent->x[1] + node->w/2.*((o/2)%2==0?1.:-1);
			node->x[2] = parent->x[2] - node->w/2.*(o>3?1.:-1);*/ 
			
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
	octant = ldexp((particles[n].x > node->x[0]) ? 0 : 1 , 0)	+ \
			ldexp((particles[n].y > node->x[1]) ? 0 : 1 , 1) + \
			ldexp((particles[n].z > node->x[2]) ? 0 : 1 , 2);
	return octant;
}

int isInside(struct cell *node){
	if (fabs(particles[node->pt].x-node->x[0]) > node->w || \
		fabs(particles[node->pt].y-node->x[1]) > node->w || \
		fabs(particles[node->pt].z-node->x[2]) > node->w) {
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
		for (int i=0; i<3; i++) {
			node->mx[i] = 0;
		}
		for (int o = 0; o < 8; o++) {
			if (node->oct[o] != NULL) {
				// Calculate the total mass and the center of mass
				for (int i=0; i<3; i++) {
					node->mx[i] = (node->mx[i]*node->m + node->oct[o]->mx[i]* \
									node->oct[o]->m) / (node->m + node->oct[o]->m);
				}
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
		root = add_particle(root, node->pt, NULL, -1);
		free(node);
		return NULL; 
	} else {
		return node;
	}
}
