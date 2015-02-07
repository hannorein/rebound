/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 		Akihiko Fujii <akihiko.fujii@nao.ac.jp>
 * @detail 	This example is identical to the shearing_sheet
 * example but uses a different algorithm for resolving individual 
 * collisions. In some cases, this might give more realistic results.
 * Particle properties resemble those found in Saturn's rings. 
 *
 * In this collision resolve method, particles are displaced if they 
 * overlap. This example also shows how to implement your own collision
 * routine. This is where one could add fragmentation, or merging of
 * particles.
 *
 * 
 * @section 	LICENSE
 * Copyright (c) 2014 Hanno Rein, Shangfei Liu, Akihiko Fujii
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
#include "collision_resolve.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "tools.h"
#include "display.h"

extern double OMEGA;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 
void collision_resolve_hardsphere_pullaway(struct collision c);

extern double opening_angle2;

void problem_init(int argc, char* argv[]){
	// Setup constants
#ifdef GRAVITY_TREE
	opening_angle2	= .5;
#endif // GRAVITY_TREE
	OMEGA 				= 0.00013143527;	// 1/s
	G 				= 6.67428e-11;		// N / (1e-5 kg)^2 m^2
	softening 			= 0.1;			// m
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
	root_nx = 2; root_ny = 2; root_nz = 1;
	nghostx = 2; nghosty = 2; nghostz = 0; 			// Use two ghost rings
	double surfacedensity 		= 840; 			// kg/m^2
	double particle_density		= 900;			// kg/m^3
	double particle_radius_min 	= 1.;			// m
	double particle_radius_max 	= 1.;			// m
	double particle_radius_slope 	= -3;	
	boxsize 			= 40;			// m
	if (argc>1){						// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}
	init_box();
	
	// Initial conditions
	printf("Toomre wavelength: %f\n",4.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G);
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	// Change collision_resolve routing from default.
	collision_resolve = collision_resolve_hardsphere_pullaway;
	// Small residual velocity to avoid particles from sinking into each other.
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear
	double total_mass = surfacedensity*boxsize_x*boxsize_y;
	double mass = 0;
	while(mass<total_mass){
		struct particle pt;
		pt.x 		= tools_uniform(-boxsize_x/2.,boxsize_x/2.);
		pt.y 		= tools_uniform(-boxsize_y/2.,boxsize_y/2.);
		pt.z 		= tools_normal(1.);					// m
		pt.vx 		= 0;
		pt.vy 		= -1.5*pt.x*OMEGA;
		pt.vz 		= 0;
		pt.ax 		= 0;
		pt.ay 		= 0;
		pt.az 		= 0;
		double radius 	= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
#ifndef COLLISIONS_NONE
		pt.r 		= radius;						// m
#endif
		double		particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
		pt.m 		= particle_mass; 	// kg
		particles_add(pt);
		mass += particle_mass;
	}
}

double coefficient_of_restitution_bridges(double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void problem_output(){

#ifdef LIBPNG
	if (output_check(2.*1e-2*2.*M_PI/OMEGA)){
		output_png("png/");
	}
#endif //LIBPNG
	if (output_check(1e-3*2.*M_PI/OMEGA)){
		output_timing();
		//output_append_velocity_dispersion("veldisp.txt");
	}
	if (output_check(2.*M_PI/OMEGA)){
		//output_ascii("position.txt");
	}
}

void problem_finish(){
}


// Function written by Akihiko Fujii
void collision_resolve_hardsphere_pullaway(struct collision c){
#ifndef COLLISIONS_NONE
	struct particle p1 = particles[c.p1];
	struct particle p2;
#ifdef MPI
	int isloc = communication_mpi_rootbox_is_local(c.ri);
	if (isloc==1){
#endif // MPI
		p2 = particles[c.p2];
#ifdef MPI
	}else{
		int root_n_per_node = root_n/mpi_num;
		int proc_id = c.ri/root_n_per_node;
		p2 = particles_recv[proc_id][c.p2];
	}
#endif // MPI
	//	if (p1.lastcollision==t || p2.lastcollision==t) return;
	struct ghostbox gb = c.gb;
	double x21  = p1.x + gb.shiftx  - p2.x; 
	double y21  = p1.y + gb.shifty  - p2.y; 
	double z21  = p1.z + gb.shiftz  - p2.z; 
	double r = sqrt(x21*x21 + y21*y21 + z21*z21);
	/* double r21 = sqrt(x21*x21 + y21*y21 + z21*z21); */
	double rp   = p1.r+p2.r;
	double oldvyouter;

	if (x21>0){
		oldvyouter = p1.vy;
	}else{
		oldvyouter = p2.vy;
	}

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
	double vy21nn = -sphi* vx21  + cphi * vy21n;

	// Coefficient of restitution
	double eps= coefficient_of_restitution_for_velocity(vx21nn);
	double dvx2 = -(1.0+eps)*vx21nn;
	double dvy2 = (r/rp-1.)*vy21nn;

	double minr = (p1.r>p2.r)?p2.r:p1.r;
	double maxr = (p1.r<p2.r)?p2.r:p1.r;
	double mindv= minr*minimum_collision_velocity;
	mindv *= 1.-(r - maxr)/minr;
	if (mindv>maxr*minimum_collision_velocity)mindv = maxr*minimum_collision_velocity;
	if (dvx2<mindv) dvx2 = mindv;

	// added
	double dxx2 = rp-r;
	double dxx2n = cphi * dxx2;
	double dxy2n = sphi * dxx2;
	double dxy2nn = ctheta * dxy2n;
	double dxz2nn = stheta * dxy2n;

	// Now we are rotating backwards
	/* double dvx2n = cphi * dvx2;	 */
	/* double dvy2n = sphi * dvx2; */

	// updated
	double dvx2n = cphi * dvx2 - sphi * dvy2;
	double dvy2n = sphi * dvx2 + cphi * dvy2;

	double dvy2nn = ctheta * dvy2n;
	double dvz2nn = stheta * dvy2n;

	// Applying the changes to the particles.
#ifdef MPI
	if (isloc==1){
#endif // MPI

		const double p1pf = p1.m/(p1.m+p2.m);
		const double p2pf = p2.m/(p1.m+p2.m);
		particles[c.p2].vx -=	p1pf*dvx2n;
		particles[c.p2].vy -=	p1pf*dvy2nn;
		particles[c.p2].vz -=	p1pf*dvz2nn;
		particles[c.p2].lastcollision = t;

		// added
		particles[c.p2].x -=	p1pf*dxx2n;
		particles[c.p2].y -=	p1pf*dxy2nn;
		particles[c.p2].z -=	p1pf*dxz2nn;
#ifdef MPI
	}
#endif // MPI
	particles[c.p1].vx +=	p2pf*dvx2n;
	particles[c.p1].vy +=	p2pf*dvy2nn;
	particles[c.p1].vz +=	p2pf*dvz2nn;

	// added
	particles[c.p1].x +=	p2pf*dxx2n; 
	particles[c.p1].y +=	p2pf*dxy2nn; 
	particles[c.p1].z +=	p2pf*dxz2nn; 

	particles[c.p1].lastcollision = t;

#endif // COLLISIONS_NONE
}

