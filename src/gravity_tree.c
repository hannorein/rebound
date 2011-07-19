#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"

double opening_angle2 = 0.25;

void calculate_forces_for_particle(int pt, struct ghostbox gb);
void calculate_forces_for_particle_from_cell(int pt, struct cell *node, struct ghostbox gb);

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
		// Summing over all particle pairs
		for (int i=0; i<N; i++){
			struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
			// Precalculated shifted position
			gb.shiftx += particles[i].x;
			gb.shifty += particles[i].y;
			gb.shiftz += particles[i].z;
			calculate_forces_for_particle(i, gb);
		}
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
	double dx = gb.shiftx - node->x;
	double dy = gb.shifty - node->y;
	double dz = gb.shiftz - node->z;
	double r2 = dx*dx + dy*dy + dz*dz;
	if (((node->w*node->w)> opening_angle2*r2) && (node->pt < 0)) {
		for (int o=0; o<8; o++) {
			if (node->oct[o] != NULL) {
				calculate_forces_for_particle_from_cell(pt, node->oct[o], gb);
			}
		}
	} else {
		if (node->pt == pt) return;
		double dx = gb.shiftx - node->mx;
		double dy = gb.shifty - node->my;
		double dz = gb.shiftz - node->mz;
		double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
		double prefact = -G/(r*r*r)*node->m;
		particles[pt].ax += prefact*dx; 
		particles[pt].ay += prefact*dy; 
		particles[pt].az += prefact*dz; 
	}
}
