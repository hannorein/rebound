/**
 * @file 	problem.c
 * @brief 	Example problem: shearing sheet.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This problem uses shearing sheet boundary
 * conditions. Particle properties resemble those found in 
 * Saturn's rings. 
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

extern double coefficient_of_restitution;
int N_border;

void problem_init(int argc, char* argv[]){
	// Setup constants
	double radius 	= 1;
	double mass	= 1;
	dt 		= 1e-1;	
	coefficient_of_restitution 	= 0.5;
	root_nx = 1; root_ny = 1; root_nz = 1;
	nghostx = 1; nghosty = 1; nghostz = 0; 	
	boxsize = 40;
	init_box();
	double N_part 	= 600;

	// Add Border Particles
	double border_spacing_x = boxsize_x/(floor(boxsize_x/radius/2.)-1.);
	printf("x: %f",border_spacing_x);
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
