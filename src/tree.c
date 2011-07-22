#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "boundaries.h"
#include "tree.h"

#if defined(GRAVITY_TREE) || defined(COLLISIONS_TREE)

struct cell** root;

int tree_get_octant_for_particle_in_cell(int pt, struct cell *node);
struct cell *tree_add_particle_to_cell(struct cell *node, int pt, struct cell *parent, int o);

void tree_add_particle_to_tree(int pt){
	if (root==NULL){
		root = calloc(root_nx*root_ny*root_nz,sizeof(struct cell*));
	}
	int root_index = get_rootbox_for_particle_int(pt);
	root[root_index] = tree_add_particle_to_cell(root[root_index],pt,NULL,0);
}

struct cell *tree_add_particle_to_cell(struct cell *node, int pt, struct cell *parent, int o){
	if (node == NULL) {
		node = calloc(1, sizeof(struct cell));
		struct particle p = particles[pt];
#ifdef GRAVITY_TREE
		node->m	 = p.m;
		node->mx = p.x;
		node->my = p.y;
		node->mz = p.z;
#endif // GRAVITY_TREE
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
	int test = -1;
	if (node == NULL) {
		return NULL;
	}
	// Non-leaf nodes	
	if (node->pt < 0) {
		for (int o = 0; o < 8; o++) {
			node->oct[o] = tree_update_cell(node->oct[o]);
		}
		// Calculate the total mass (, quadrupole tensor) and center of mass, and check
		// if the node needs derefinement after updating the tree.
		node->pt = 0;
#ifdef GRAVITY_TREE
		node->m	 = 0;
		node->mx = 0;
		node->my = 0;
		node->mz = 0;
	#ifdef QUADRUPOLE
		node->mxx = 0;
		node->mxy = 0;
		node->mxz = 0;
		node->myy = 0;
		node->myz = 0;
		node->mzz = 0;
	#endif // QUADRUPOLE
#endif // GRAVITY_TREE
		for (int o = 0; o < 8; o++) {
			if (node->oct[o] != NULL) {
				// Calculate the total mass and the center of mass
#ifdef GRAVITY_TREE
				double m_o = node->oct[o]->m;
				node->mx = (node->mx*node->m + node->oct[o]->mx*m_o) / (node->m + m_o);
				node->my = (node->my*node->m + node->oct[o]->my*m_o) / (node->m + m_o);
				node->mz = (node->mz*node->m + node->oct[o]->mz*m_o) / (node->m + m_o);
				node->m += m_o;
	#ifdef QUADRUPOLE
				// Ref: Hernquist, L., 1987, APJS
				double qx  = node->oct[o]->mx - node->mx;
				double qy  = node->oct[o]->my - node->my;
				double qz  = node->oct[o]->mz - node->mz;
				double qr2 = qx*qx + qy*qy + qz*qz;
				node->mxx += node->oct[o]->mxx + m_o*(3*qx*qx - qr2);
				node->mxy += node->oct[o]->mxy + m_o*3*qx*qy;
				node->mxz += node->oct[o]->mxz + m_o*3*qx*qz;
				node->myy += node->oct[o]->myy + m_o*(3*qy*qy - qr2);
				node->myz += node->oct[o]->myz + m_o*3*qy*qz;
				node->mzz += -node->mxx -node->myy;
	#endif
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
#ifdef GRAVITY_TREE
		struct particle p = particles[node->pt];
		node->mx = p.x;
		node->my = p.y;
		node->mz = p.z;
#endif // GRAVITY_TREE
		particles[node->pt].c = node;
		return node;
	}
}

void tree_update(){
	if (root==NULL){
		root = calloc(root_nx*root_ny*root_nz,sizeof(struct cell*));
	}
	for(int i=0;i<root_nx;i++){
	for(int j=0;j<root_ny;j++){
	for(int k=0;k<root_nz;k++){
		int index = (k*root_ny+j)*root_nx+i;
		root[index] = tree_update_cell(root[index]);
	}
	}
	}
}


#endif
