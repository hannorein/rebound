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
void calculate_forces_for_particle_from_cell(int pt, struct cell const *node, struct ghostbox const gb);

void calculate_forces(){
#pragma omp parallel for
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
#pragma omp parallel for
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
		struct cell* node = tree_root[i];
		if (node!=NULL){
			calculate_forces_for_particle_from_cell(pt, node, gb);
		}
	}
}

void calculate_forces_for_particle_from_cell(int pt, struct cell const *node, struct ghostbox const gb) {
	double dx = gb.shiftx - node->mx;
	double dy = gb.shifty - node->my;
	double dz = gb.shiftz - node->mz;
	double r2 = dx*dx + dy*dy + dz*dz;
	if (((node->w*node->w)> opening_angle2*r2) && (node->pt < 0)) {
		for (int o=0; o<8; o++) {
			if (node->oct[o] != NULL) {
				calculate_forces_for_particle_from_cell(pt, node->oct[o], gb);
			}
		}
	} else {
		if (node->pt == pt) return;
		double r = sqrt(r2 + softening*softening);
		double prefact = -G/(r*r*r)*node->m;
#ifdef QUADRUPOLE
		if (node->pt < 0) {
			double qprefact = G/(r*r*r*r*r);
			particles[pt].ax += qprefact*(dx*node->mxx + dy*node->mxy + dz*node->mxz); 
			particles[pt].ay += qprefact*(dx*node->mxy + dy*node->myy + dz*node->myz); 
			particles[pt].az += qprefact*(dx*node->mxz + dy*node->myz + dz*node->mzz); 
			double mrr 	= dx*dx*node->mxx 	+ dy*dy*node->myy 	+ dz*dz*node->mzz
					+ 2.*dx*dy*node->mxy 	+ 2.*dx*dz*node->mxz 	+ 2.*dy*dz*node->myz; 
			qprefact *= -5.0/(2.0*r*r)*mrr;
			particles[pt].ax += qprefact*dx; 
			particles[pt].ay += qprefact*dy; 
			particles[pt].az += qprefact*dz; 
		}
#endif
		particles[pt].ax += prefact*dx; 
		particles[pt].ay += prefact*dy; 
		particles[pt].az += prefact*dz; 
	}
}
