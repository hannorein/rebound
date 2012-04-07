/**
 * @file 	problem.c
 * @brief 	Example problem: Granular dynamics.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses two boundary layers made of 
 * particles to simulate shearing walls. These walls are heating
 * up the particles, create a dense and cool layer in the middle.
 * This example shows how to use REBOUND for granular dynamics.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "collision_resolve.h"

void collision_resolve_hardsphere_withborder(struct collision c);
extern double coefficient_of_restitution;
int N_border;

void problem_init(int argc, char* argv[]){
	// Setup constants
	double radius 	= 1;
	double mass	= 1;
	dt 		= 1e-1;	
	// Override default collision handling to account for border particles
	coefficient_of_restitution 	= 0.15;
	collision_resolve = collision_resolve_hardsphere_withborder;
	root_nx = 1; root_ny = 1; root_nz = 4;
	nghostx = 1; nghosty = 1; nghostz = 0; 	
	boxsize = 20;
	init_box();
	double N_part 	= 0.00937*boxsize_x*boxsize_y*boxsize_z;

	// Add Border Particles
	double border_spacing_x = boxsize_x/(floor(boxsize_x/radius/2.)-1.);
	double border_spacing_y = boxsize_y/(floor(boxsize_y/radius/2.)-1.);
	struct particle pt;
	pt.vx 		= 0;
	pt.vz 		= 0;
	pt.ax 		= 0;
	pt.ay 		= 0;
	pt.az 		= 0;
	pt.r 		= radius;
	pt.m 		= mass;
	for(double x = -boxsize_x/2.; x<boxsize_x/2.-border_spacing_x/2.;x+=border_spacing_x){
		for(double y = -boxsize_y/2.; y<boxsize_y/2.-border_spacing_y/2.;y+=border_spacing_y){
			pt.x 		= x;
			pt.y 		= y;
			
			// Add particle to bottom
			pt.z 		= -boxsize_z/2.+radius;
			pt.vy 		= 1;
			particles_add(pt);

			// Add particle to top
			pt.z 		= boxsize_z/2.-radius;
			pt.vy 		= -1;
			particles_add(pt);
			N_tree_fixed++;
		}
	}

	N_border = N;
#ifdef TREE
	N_tree_fixed = N_border;
#endif // TREE
	
	// Add real particles
	while(N-N_border<N_part){
		struct particle pt;
		pt.x 		= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
		pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
		pt.z 		= 0.758*tools_uniform(-boxsize_z/2.,boxsize_z/2.);
		pt.vx 		= tools_normal(0.001);
		pt.vy 		= tools_normal(0.001);
		pt.vz 		= tools_normal(0.001);
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		pt.r 		= radius;						// m
		pt.m 		= 1;
		particles_add(pt);
	}
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(10*dt)){
		output_timing();
	}
}

void problem_finish(){
}

void collision_resolve_hardsphere_withborder(struct collision c){
	struct particle p1 = particles[c.p1];
	struct particle p2 = particles[c.p2];
	struct ghostbox gb = c.gb;
	double m21  = p1.m  /  p2.m; 
	double x21  = p1.x + gb.shiftx  - p2.x; 
	double y21  = p1.y + gb.shifty  - p2.y; 
	double z21  = p1.z + gb.shiftz  - p2.z; 
	double rp   = p1.r+p2.r;
	if (rp*rp < x21*x21 + y21*y21 + z21*z21) return;
	double vx21 = p1.vx + gb.shiftvx - p2.vx; 
	double vy21 = p1.vy + gb.shiftvy - p2.vy; 
	double vz21 = p1.vz + gb.shiftvz - p2.vz; 
	if (vx21*x21 + vy21*y21 + vz21*z21 >0) return; // not approaching
	// Bring the to balls in the xy plane.
	// NOTE: this could probabely be an atan (which is faster than atan2)
	double theta = atan2(z21,y21);
	double stheta = sin(theta);
	double ctheta = cos(theta);
	double vy21n = ctheta * vy21 + stheta * vz21;	
	double y21n = ctheta * y21 + stheta * z21;	
	
	// Bring the two balls onto the positive x axis.
	double phi = atan2(y21n,x21);
	double cphi = cos(phi);
	double sphi = sin(phi);
	double vx21nn = cphi * vx21  + sphi * vy21n;		

	// Coefficient of restitution
	double eps= coefficient_of_restitution_for_velocity(vx21nn);
	double dvx2 = -(1.0+eps)*vx21nn/(1.0+m21) ;
	if (dvx2<minimum_collision_velocity){
		dvx2 = minimum_collision_velocity;
	}

	// Now we are rotating backwards
	double dvx2n = cphi * dvx2;		
	double dvy2n = sphi * dvx2;		
	double dvy2nn = ctheta * dvy2n;	
	double dvz2nn = stheta * dvy2n;	

	// Applying the changes to the particles.
	if (c.p2>N_border){
		particles[c.p2].vx -=	m21*dvx2n;
		particles[c.p2].vy -=	m21*dvy2nn;
		particles[c.p2].vz -=	m21*dvz2nn;
		particles[c.p2].lastcollision = t;
	}
	if (c.p1>N_border){
		particles[c.p1].vx +=	dvx2n; 
		particles[c.p1].vy +=	dvy2nn; 
		particles[c.p1].vz +=	dvz2nn; 
		particles[c.p1].lastcollision = t;
	}
}
