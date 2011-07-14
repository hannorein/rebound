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
/*	int	level;
	level = 0;*/
	root = NULL;
	for (int i=0;i<N;i++){
//		root = add_leaf(root, particles[i], level);
		root = add_leaf(root, i, NULL, -1);
	}
}

struct cell *add_leaf(struct cell *node, int n, struct cell *parent, int o){
	if (node == NULL) {
		node	= malloc(sizeof(struct cell));
//		node.m	= particles[n].m;
//		node.mx	= particles[n].x;
//		node.my	= particles[n].y;
//		node.mz	= particles[n].z;
		if (parent == NULL){
			for (int i=0; i<3; i++) node->x[i] = 0;
		/*	node->x	= 0;
			node->y	= 0;
			node->z	= 0; */
			node->w = boxsize;
		} else {
			node->w = parent->w/2.;
//			node->x	= node->w * (floor(node->mx / node->w) + 0.5);
//			node->y	= node->w * (floor(node->my / node->w) + 0.5);
//			node->z	= node->w * (floor(node->mz / node->w) + 0.5);
			for (int i=0; i<3; i++) {
				node->x[i] = parent->x[i] + node->w/2.*((o>>i)%2==0?1.:-1);
			}
		/*	node->x[0] = parent->x[0] + node->w/2.*(o%2==0?1.:-1);
			node->x[1] = parent->x[1] + node->w/2.*((o/2)%2==0?1.:-1);
			node->x[2] = parent->x[2] - node->w/2.*(o>3?1.:-1);*/ 
			
		}
		node->pt	= n;	//?
//		node->lv	= l; 
//		node->num = 1;
//		particles[n].leaf = node; 
		for (int i=0;i<8;i++) node->oct[i] = NULL;
		//printf("%d%",node->w);
	} else if (isLeaf(node->oct) == 1){
//		node->num++;
//		node.m	+= particles[n].m;
//		node.mx	=  (node.mx*node.m + particles[n].m*particles[n].x) / node.m;
//		node.my	=  (node.my*node.m + particles[n].m*particles[n].y) / node.m;
//		node.mz	=  (node.mz*node.m + particles[n].m*particles[n].z) / node.m;
//		if (node->num == 2){
		int o;
		o = get_octant(node->pt, node);
		node->oct[o] = add_leaf(node->oct[o], node->pt, node, o); 
		o = get_octant(n, node);
		node->oct[o] = add_leaf(node->oct[o], n, node, o);
//			particles[node->pt].leaf = NULL;
		}
	else {
		int o = get_octant(n, node);
		node->oct[o] = add_leaf(node->oct[o], n, node, o);
	}
	return node;
}

int get_octant(int n, struct cell *node){
	int octant = 0;
/*	for (int i=0; i<3; i++) {
		octant += ldexp((particles[n].x[i] > node->x[i]) ? 0 : 1 , 0);
	}*/
	octant = ldexp((particles[n].x > node->x[0]) ? 0 : 1 , 0)	+ \
			ldexp((particles[n].y > node->x[1]) ? 0 : 1 , 1) + \
			ldexp((particles[n].z > node->x[2]) ? 0 : 1 , 2);
	return octant;
}

int isLeaf(struct cell **octant){
	for (int o=0; o<8; o++) {
	//	if (*(octant+o) != NULL) return 0;
		if (octant[o] != NULL) return 0;
	}
	return 1;
}

int need_derefine(struct cell *node){
	int leaves = 0;
	int oct = -1;
	for (int o=0; o<8; o++) {
		if (node->oct[o]!=NULL){
			if ((isLeaf(node->oct[o]->oct) == 0) || leaves > 1){
				return 0;
			}else{
				leaves++;
				oct = o;
			}
		}
	}
	return oct;
}

int isInside(struct cell *node){
/*	for (int i=0; i<3; i++) {
		if (fabs(particles[node->pt].x[i]-node->x[i]) > node->w) {
			return 0;
		}
	}*/
	if (fabs(particles[node->pt].x-node->x[0]) > node->w || \
		fabs(particles[node->pt].y-node->x[1]) > node->w || \
		fabs(particles[node->pt].z-node->x[2]) > node->w) {
		return 0;
	}
	return 1;
}
	
/*
void cell *derefine(struct cell *node) {
	int o = get_octant(node->pt, node);
	for (int i = 0; i<3; i++) {
		node->x[i] = node->oct[o]->x[i];
	}
	free(node->oct[o]);
}*/

struct cell *update_tree(struct cell *node){
	//printf("Updating the tree\n");
/*	for (int i=0;i<N;i++) {
		if (fabs(particles[i].x-particles[i].tree->x) > particles[i].tree->w || \
			fabs(particles[i].y-particles[i].tree->y) > particles[i].tree->w || \
			fabs(particles[i].z-particles[i].tree->z) > particles[i].tree->w ) {
				particles[i].tree->
				free();	
				add_leaf(root, i, NULL, -1);
		} 
	} */
	int test;
	if (node == NULL) {
		return NULL;
	} else	if (isLeaf(node->oct) == 0) {
//		int leaves = 0;
		for (int o = 0; o < 8; o++) {
			node->oct[o] = update_tree(node->oct[o]);
//			if (node->oct[o] == NULL) leaves++;
		}
		// If no leaf left
/*		if (leaves == 0) {
			free(node);
			return NULL;
		} */
		// Check if needs derefinement.
		if (isLeaf(node->oct) == 1) {
			free(node);
			return NULL;
		} else if ((test = need_derefine(node)) != 0) {
			node->pt = node->oct[test]->pt;
			free(node->oct[test]);
			node->oct[test]=NULL;
			return node;
		}
		return node;
	} else if (isInside(node) == 0) { 
		root = add_leaf(root, node->pt, NULL, -1);
		free(node);
		return NULL; 
	} else {
		return node;
	}
//	return node;
}
