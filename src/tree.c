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
	int	level;
	level = 0;
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
			node->x	= 0;
			node->y	= 0;
			node->z	= 0;
			node->w = boxsize;
		} else {
			node->w = parent->w/2.;
//			node->x	= node->w * (floor(node->mx / node->w) + 0.5);
//			node->y	= node->w * (floor(node->my / node->w) + 0.5);
//			node->z	= node->w * (floor(node->mz / node->w) + 0.5);
			node->x = parent->x + node->w/2.*(o%2==0?1.:-1);
			node->y = parent->y + node->w/2.*((o/2)%2==0?1.:-1);
			node->z = parent->z - node->w/2.*(o>3?1.:-1);
			
		}
		node->pt	= n;	//?
//		node->lv	= l; 
		node->num = 1;
		for (int i=0;i<8;i++) node->oct[i] = NULL;
		//printf("%d%",node->w);
	} else {
		node->num++;
//		node.m	+= particles[n].m;
//		node.mx	=  (node.mx*node.m + particles[n].m*particles[n].x) / node.m;
//		node.my	=  (node.my*node.m + particles[n].m*particles[n].y) / node.m;
//		node.mz	=  (node.mz*node.m + particles[n].m*particles[n].z) / node.m;
		if (node->num == 2){
			int o;
			o = get_octant(node->pt, node);
			node->oct[o] = add_leaf(node->oct[o], node->pt, node, o); 
			o = get_octant(n, node);
			node->oct[o] = add_leaf(node->oct[o], n, node, o);
		}
		else {
			int o = get_octant(n, node);
			node->oct[o] = add_leaf(node->oct[o], n, node, o);
		}
	}
	return node;
}

int get_octant(int n, struct cell *node){
	int octant;
	octant = ldexp((particles[n].x>node->x) ? 0 : 1 , 0)	+ \
			ldexp((particles[n].y>node->y) ? 0 : 1 , 1) + \
			ldexp((particles[n].z>node->z) ? 0 : 1 , 2);
	return octant;
}
