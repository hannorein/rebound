#ifdef GRAVITY_TREE
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

void add_particle(struct cell* node, struct particle* pt);


void tree_init(){						// Initializes the tree and adds all particles.
	printf("Initializing the tree\n");
	root = calloc(1,sizeof(struct cell));
	root->w = boxsize;
	for (int i=0;i<N;i++){
		add_particle(root, &(particles[i]));
	}
}

int get_daughter_index(struct cell *node, struct particle* pt){	// Get's the octant index for the particle pt in cell node.
	int octant = 0;
	if (pt->x>node->x[0]) octant+=1;
	if (pt->y>node->x[1]) octant+=2;
	if (pt->z>node->x[2]) octant+=4;
	return octant;
}

void add_particle_to_daughter(struct cell* node, struct particle* pt){ // Puts the particle pt in the correct daughter cell of node
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

void add_particle(struct cell* node, struct particle* pt){	// Tries to put the particle pt in the cell node
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
		if (fabs(p->x-node->x[0]) > node->w/2. || 
			fabs(p->y-node->x[1]) > node->w/2. || 
			fabs(p->z-node->x[2]) > node->w/2.) {
			node->pt=NULL;				// Removes particle from tree. Cleanup is done later.
			add_particle(root,p);			// Readds particle to tree.
		}
	}
}
void tree_update(struct cell *node){
	if (node->oct!=NULL){					// Can only simplify daughters if there are any.
		int daughters = 0;
		int nodelast = 0;
		for (int i=0;i<8;i++){
			if (node->oct[i]!=NULL){
				tree_update(node->oct[i]);  	// Recursive cleanup
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
#endif //GRAVITY_TREE
