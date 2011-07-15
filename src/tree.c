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

void add_particle_to_cell(struct cell* node, struct particle* pt);
void add_particle_to_tree(struct particle* pt);


void tree_init(){						// Initializes the tree and adds all particles.
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
		root[index] = calloc(1,sizeof(struct cell));
		root[index]->w = boxsize_min;
		root[index]->x = -boxsize_x/2.+boxsize_min*(0.5+(double)i);
		root[index]->y = -boxsize_y/2.+boxsize_min*(0.5+(double)j);
		root[index]->z = -boxsize_z/2.+boxsize_min*(0.5+(double)k);
	}
	}
	}

	for (int i=0;i<N;i++){
		add_particle_to_tree(&(particles[i]));
	}
}

int get_daughter_index(struct cell *node, struct particle* pt){	// Get's the octant index for the particle pt in cell node.
	int octant = 0;
	if (pt->x>node->x) octant+=1;
	if (pt->y>node->y) octant+=2;
	if (pt->z>node->z) octant+=4;
	return octant;
}

void add_particle_to_daughter(struct cell* node, struct particle* pt){ // Puts the particle pt in the correct daughter cell of node
	int o = get_daughter_index(node,pt);
	if (node->oct[o]==NULL){				
		node->oct[o]=calloc(1,sizeof(struct cell));
		node->oct[o]->w = node->w/2.; 
		node->oct[o]->x = node->x + node->w/4.*((o>>0)%2==0?-1.:+1.);
		node->oct[o]->y = node->y + node->w/4.*((o>>1)%2==0?-1.:+1.);
		node->oct[o]->z = node->z + node->w/4.*((o>>2)%2==0?-1.:+1.);
	}
	add_particle_to_cell(node->oct[o],pt);

}

void add_particle_to_tree(struct particle* pt){				// Tries to put the particle pt in tree
	int i = ((int)floor((pt->x + boxsize_x/2.)/boxsize_min));
	int j = ((int)floor((pt->y + boxsize_y/2.)/boxsize_min));
	int k = ((int)floor((pt->z + boxsize_z/2.)/boxsize_min));
	int index = (k*root_ny+j)*root_nx+i;
	add_particle_to_cell(root[index],pt);
}

void add_particle_to_cell(struct cell* node, struct particle* pt){	// Tries to put the particle pt in the cell node
	if (node->oct==NULL){
		if(node->pt==NULL){				// Leaf with no particle in it. Put pt here.
			node->pt = pt;
			pt->c = node;
			return;
		}else{						// Leaf with one particle in it. 
			node->oct = calloc(8,sizeof(struct cell*));
			add_particle_to_daughter(node,node->pt);// Then re-add particle.
			node->pt=NULL;				// Need to remove particle in cell.
		}
	}
	add_particle_to_daughter(node,pt);			// Not a leaf, put particle in daughter cell.
}
	
void tree_check_moved_particles(){				// Check each particle if it's in it's cell.
	for (int i=0;i<N;i++){
		struct particle* p = &(particles[i]);
		struct cell* node = p->c;
		if (fabs(p->x-node->x) > node->w/2. || 
			fabs(p->y-node->y) > node->w/2. || 
			fabs(p->z-node->z) > node->w/2.) {
			node->pt=NULL;				// Removes particle from tree. Cleanup is done later.
			add_particle_to_tree(p);			// Readds particle to tree.
		}
	}
}
void tree_update(){
	if (root==NULL){
		tree_init();
	} else {
		check_boundaries();
		tree_check_moved_particles();
		for(int i=0;i<root_nx;i++){
		for(int j=0;j<root_ny;j++){
		for(int k=0;k<root_nz;k++){
			int index = (k*root_ny+j)*root_nx+i;
			tree_update_cell(root[index]);
		}
		}
		}
	}
}

void tree_update_cell(struct cell *node){
	if (node->oct!=NULL){					// Can only simplify daughters if there are any.
		int daughters = 0;
		int nodelast = 0;
		for (int i=0;i<8;i++){
			if (node->oct[i]!=NULL){
				tree_update_cell(node->oct[i]);  	// Recursive cleanup
				if (node->oct[i]->pt==NULL){
					if (node->oct[i]->oct!=NULL){
						daughters=2;	// cannot simplify further
					}
				}else{
					daughters++;		// keep track of number of daughters with one particle inside
					nodelast=i;
				}
			}
		}
		if (daughters==0){				// No daughters. Delete oct.
			free(node->oct);
			node->oct=NULL;
		}
		if (daughters==1){				// Exactly one daughter with one particle. Can simplify.
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
