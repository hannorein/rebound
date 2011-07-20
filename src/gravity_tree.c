#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "particle.h"
#include "main.h"
#include "tree.h"
#include "boundaries.h"
#ifndef QUADRUPOLE
	#define QUADRUPOLE
#endif

double opening_angle2 = 0.25;

void calculate_forces_for_particle(int pt, struct ghostbox gb);
void calculate_forces_for_particle_from_cell(int pt, struct cell *node, struct ghostbox gb);
void calculate_quadrupole_tensor(double *mqp, struct cell *node, double mcx, double mcy, double mcz); 

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
		particles[pt].ax += prefact*dx; 
		particles[pt].ay += prefact*dy; 
		particles[pt].az += prefact*dz; 
#ifdef QUADRUPOLE
		if (node->pt < 0) {
			double mq[3][3];
			for (int i=0; i<9; i++) {
				mq[i/3][i%3] = 0;
			}
			calculate_quadrupole_tensor(&mq[0][0], node, node->mx, node->my, node->mz);
			double m4 = dx*dx*mq[0][0] + \
						dx*dy*mq[0][1] + \
						dx*dz*mq[0][2] + \
						dy*dx*mq[1][0] + \
						dy*dy*mq[1][1] + \
						dy*dz*mq[1][2] + \
						dz*dx*mq[2][0] + \
						dz*dy*mq[2][1] + \
						dz*dz*mq[2][2]; 
			double qprefact = -G/(2*r*r*r*r*r)*m4;
			particles[pt].ax += qprefact*dx; 
			particles[pt].ay += qprefact*dy; 
			particles[pt].az += qprefact*dz; 
		}
#endif
	}
}

#ifdef QUADRUPOLE
void calculate_quadrupole_tensor(double *mqp, struct cell *node, double mcx, double mcy, double mcz) {
	for (int o=0; o<8; o++) {
		if (node->oct[o] != NULL) {
			if (node->oct[o]->pt >= 0) {
				double pm	= node->oct[o]->m;
				double px	= node->oct[o]->mx - mcx;
				double py	= node->oct[o]->my - mcy;
				double pz	= node->oct[o]->mz - mcz;
				double pr2	= px*px + py*py + pz*pz;
				*(mqp)   += pm * (3*px*px - pr2);
				*(mqp+1)   += pm * 3*px*py; 
				*(mqp+2)   += pm * 3*px*pz; 
				*(mqp+3)   += pm * 3*py*px; 
				*(mqp+4)   += pm * (3*py*py - pr2);
				*(mqp+5)   += pm * 3*py*pz; 
				*(mqp+6)   += pm * 3*pz*px; 
				*(mqp+7)   += pm * 3*pz*py; 
				*(mqp+8)   += pm * (3*pz*pz - pr2);
			/*	mq[0][0]   += pm * (3*px*px - pr2);
				mq[0][1]   += pm * 3*px*py; 
				mq[0][2]   += pm * 3*px*pz; 
				mq[1][0]   += pm * 3*py*px; 
				mq[1][1]   += pm * (3*py*py - pr2);
				mq[1][2]   += pm * 3*py*pz; 
				mq[2][0]   += pm * 3*pz*px; 
				mq[2][1]   += pm * 3*pz*py; 
				mq[2][2]   += pm * (3*pz*pz - pr2); */
			} else {
				calculate_quadrupole_tensor(mqp, node->oct[o], mcx, mcy, mcz);
			}
		}
	}
}
#endif
