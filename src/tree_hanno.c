#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "tree_hanno.h"

struct cell* root;

void add_particle(struct cell* node, struct particle* pt);


void init_tree(){
	printf("Initializing the tree\n");
	root = calloc(1,sizeof(struct cell));
	root->w = boxsize;
	for (int i=0;i<N;i++){
		add_particle(root, &(particles[i]));
	}
}

int get_daughter_index(struct cell *node, struct particle* pt){
	int octant = 0;
	if (pt->x>node->x[0]) octant+=1;
	if (pt->y>node->x[1]) octant+=2;
	if (pt->z>node->x[2]) octant+=4;
	return octant;
}

void add_particle_to_daughter(struct cell* node, struct particle* pt){
	int o = get_daughter_index(node,pt);
	if (node->oct[o]==NULL){
		node->oct[o]=calloc(1,sizeof(struct cell));
		node->oct[o]->w = node->w/2.; 
		for (int i=0; i<3; i++) {
			node->oct[o]->x[i] = node->x[i] + node->w/4.*((o>>i)%2==0?-1.:+1.);
		}
	}
	add_particle(node->oct[o],pt);

}

void add_particle(struct cell* node, struct particle* pt){
	if (node->oct==NULL){
		if(node->pt==NULL){
			// Leaf with no particle in it
			node->pt = pt;
			pt->c = node;
			return;
		}else{
			// Leaf with one particle in it
			node->oct = calloc(8,sizeof(struct cell*));
			add_particle_to_daughter(node,node->pt);
			node->pt=NULL;
		}
	}
	// Not a leaf, put particle in daughter
	add_particle_to_daughter(node,pt);
}

int isInside(struct cell *node){
	if (fabs(node->pt->x-node->x[0]) > node->w/2. || 
		fabs(node->pt->y-node->x[1]) > node->w/2. || 
		fabs(node->pt->z-node->x[2]) > node->w/2.) {
		return 0;
	}
	return 1;
}
	
void update_tree1(){
	for (int i=0;i<N;i++){
		if (!isInside(particles[i].c)){
			particles[i].c->pt=NULL;
			add_particle(root,&(particles[i]));
		}
	}
}
void update_tree(struct cell *node){
	if (node->oct!=NULL){
		int daughters = 0;
		int nodelast = 0;
		for (int i=0;i<8;i++){
			if (node->oct[i]!=NULL){
				update_tree(node->oct[i]);
				if (node->oct[i]->pt==NULL){
					if (node->oct[i]->oct!=NULL){
						daughters=2;
					}
				}else{
					daughters++;
					nodelast=i;
				}
			}
		}

		if (daughters==0){
			free(node->oct);
			node->oct=NULL;
		}
		if (daughters==1){
			if(node->oct[nodelast]->pt!=NULL){
				node->pt = node->oct[nodelast]->pt;
				node->pt->c = node;
				free(node->oct[nodelast]);
				free(node->oct);
				node->oct=NULL;
			}
		}
	}
}
