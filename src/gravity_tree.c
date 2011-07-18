#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"
#ifndef THETA
	#define THETA 0.2
#endif

void calculate_forces(){
	tree_update();
	for (int i=0; i<N; i++){
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
	}
	// Summing over all Ghost Boxes
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
	for (int gbz=-nghostz; gbz<=nghostz; gbz++){
		struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
		// Summing over all particle pairs
		for (int i=0; i<N; i++){
//		for (int j=0; j<N; j++){
//			if (i==j) continue;
			calculate_forces_for_particle(i, gb);
//			double dx = particles[i].x - (particles[j].x + gb.shiftx);
//			double dy = particles[i].y - (particles[j].y + gb.shifty);
//			double dz = particles[i].z - (particles[j].z + gb.shiftz);
//			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
//			double prefact = -G/(r*r*r)*particles[j].m;
//			particles[i].ax += prefact*dx; 
//			particles[i].ay += prefact*dy; 
//			particles[i].az += prefact*dz; 
		}
//		}
	}
	}
	}
}

void calculate_forces_for_particle(int pt, struct ghostbox gb) {
	int root_n = root_nx*root_ny*root_nz;
	for(int i=0;i<root_n;i++){
		struct cell* node = root[i];
		if (node!=NULL){
			calculate_forces_for_particle_from_cell(pt, node, gb);
		}
	}
}

void calculate_forces_for_particle_from_cell(int pt, struct cell *node, struct ghostbox gb) {
	double dx = particles[pt].x - (node->x + gb.shiftx);
	double dy = particles[pt].y - (node->y + gb.shifty);
	double dz = particles[pt].z - (node->z + gb.shiftz);
	double r2 = dx*dx + dy*dy + dz*dz;
	if (((node->w*node->w)/r2 > THETA*THETA) && (node->pt < 0)) {
		for (int o=0; o<8; o++) {
			if (node->oct[o] != NULL) {
				calculate_forces_for_particle_from_cell(pt, node->oct[o], gb);
			}
		}
	} else {
		if (node->pt == pt) return;
		double dx = particles[pt].x - (node->mx + gb.shiftx);
		double dy = particles[pt].y - (node->my + gb.shifty);
		double dz = particles[pt].z - (node->mz + gb.shiftz);
		double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
		double prefact = -G/(r*r*r)*node->m;
		particles[pt].ax += prefact*dx; 
		particles[pt].ay += prefact*dy; 
		particles[pt].az += prefact*dz; 
	}
}
