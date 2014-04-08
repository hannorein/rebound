/**
 * @file 	problem.c
 * @brief 	Test problem calculating various parts of the viscosity in a planetary ring.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements a test problem similar to simulations performed by
 * Daisaka et al 2001. Translational, collisional and gravitational parts of the viscosity
 * are time averaged and outputted at the end of the simulation. See Daisaka et al for
 * parameter definitions.
 * 
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
#include "tree.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
extern double opening_angle2;
const double particle_radius	= 1;	// Physical particle radius
double r_hstar;				// Dimensionless particle radius
double lambda_crit;
double sx2;				// Number of critical wavelengths in x direction (*0.5)
double sy2;				// Number of critical wavelengths in y direction (*0.5)

void problem_init(int argc, char* argv[]){
	// Check for command line arguments
	if (argc<3){
		printf("Error. Please specify two command line arguments: r_h^* and s_x\n");
		exit(1);
	}

	// Calculating parameters. See Daisaka et al for definitions.
	r_hstar = atof(argv[1]);
	double r_h 		= r_hstar*2.*particle_radius;
	double particle_mass 	= r_h*r_h*r_h*3./2.;
	double tau 		= 0.5;
	double surface_density 	= tau/(M_PI*particle_radius*particle_radius)*particle_mass;
	lambda_crit	 	= 4.*M_PI*M_PI*G*surface_density/OMEGA/OMEGA;
	sx2 			= 3.*lambda_crit/2.;
	sy2 			= 3.*lambda_crit/2.;
	boxsize 		= lambda_crit*atof(argv[2]);
	int _N 			= (int)round(tau*boxsize*boxsize/(M_PI*particle_radius*particle_radius));
	dt 			= 2e-3*2.*M_PI;
	tmax			= 40.*2.*M_PI;
	softening 		= 0.1*particle_radius;
	opening_angle2		= 0.9;
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 2.*dt*particle_radius;
	nghostx = 3; nghosty = 3; nghostz = 0;
	init_box();

	// Initialize particles
	while(N<_N){
		struct particle p;
		p.x 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z 	= 3.*particle_radius*((double)rand()/(double)RAND_MAX-0.5);
		p.vx 	= 0;
		p.vy 	= -1.5*p.x*OMEGA;
		p.vz 	= 0;
		p.ax 	= 0;
		p.ay 	= 0;
		p.az 	= 0;
		p.m 	= particle_mass;
		p.r 	= particle_radius;
		particles_add(p);
	}
}

// This function calculates the translational part of the viscosity
double calculate_viscosity_trans(){
	double nu = 0;
	double mtot = 0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		mtot	+= p.m;
		double vyr = p.vy + 1.5 * OMEGA * p.x;
		nu	+= p.m*p.vx*vyr;
	}
	return 2./(3.*OMEGA*mtot)*nu;
}

extern double collisions_plog;
double collisions_plog_last =0;

// This function calculates the collisional part of the viscosity
double calculate_viscosity_coll(){
	if (t-collisions_plog_last==0 || collisions_plog==0){
		return 0;
	}
	double mtot = 0;
	for (int i=0;i<N;i++){
		mtot	+= particles[i].m;
	}
	double tmp_dt 		= (t-collisions_plog_last);
	double tmp_plog 	= collisions_plog;
	return 2./(3.*OMEGA*mtot*tmp_dt)*tmp_plog;
}
		
// These functions calculate the gravitational part of the viscosity using the tree
double calculate_viscosity_grav_for_particle_from_cell(int i, struct ghostbox gb,struct cell* node){
	double nu = 0;
	double dx = gb.shiftx - node->mx;
	double dy = gb.shifty - node->my;
	double dz = gb.shiftz - node->mz;
	if (node->pt<0){
		if (fabs(dx)-node->w/2. < sx2 && fabs(dy)-node->w/2. < sy2){
			for(int o=0;o<8;o++){
				if (node->oct[o]!=NULL){
					nu += calculate_viscosity_grav_for_particle_from_cell(i,gb,node->oct[o]);
				}
			}
		} 
	}else{
		if (fabs(dx) < sx2 && fabs(dy) < sy2 && dx > 0 && i!=node->pt){
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			nu += particles[i].m*particles[node->pt].m*dy*dx/(r*r*r);
		}
	}
	return nu;
}

double calculate_viscosity_grav_for_particle(int pt, struct ghostbox gb) {
	double nu =0;
	for(int i=0;i<root_n;i++){
		struct cell* node = tree_root[i];
		if (node!=NULL){
			nu += calculate_viscosity_grav_for_particle_from_cell(pt, gb, node);
		}
	}
	return nu;
}

double calculate_viscosity_grav(){
	double nu  	= 0;
	double mtot 	= 0;
	for (int i=0;i<N;i++){
		mtot	+= particles[i].m;
	}

	// Summing over all Ghost Boxes
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
		struct ghostbox gb_cache = boundaries_get_ghostbox(gbx,gby,0);
		// Summing over all particle pairs
		for (int i=0; i<N; i++){
			struct ghostbox gb = gb_cache;
			// Precalculated shifted position
			gb.shiftx += particles[i].x;
			gb.shifty += particles[i].y;
			gb.shiftz += particles[i].z;
			nu += calculate_viscosity_grav_for_particle(i, gb);
		}
	}}
	return 2./(3.*OMEGA*mtot)*nu * (-G);
}


double nutrans_tot=0, nucoll_tot=0, nugrav_tot=0; int nuN=0;

int reset = 1;
void problem_output(){
	output_timing();
	if (t>tmax/2.){
      		if(reset==1){
			collisions_plog_last 	= t;
			collisions_plog		= 0;
			nuN 			= 0;
			nutrans_tot 		= 0;
			nucoll_tot  		= 0; 
			nugrav_tot  		= 0;
			reset 			= 0;
		}else{
			// Calculate viscosity
			nutrans_tot 	+= calculate_viscosity_trans();
			nugrav_tot 	+= calculate_viscosity_grav();
			nucoll_tot	 = calculate_viscosity_coll();
			nuN++;
		}
	}
}


void problem_finish(){
	// Output viscosity before exiting.
	FILE* ofv = fopen("viscosity.txt","a"); 
	double nu_trans = nutrans_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN;
	double nu_coll 	= nucoll_tot/(particle_radius*particle_radius*OMEGA);
	double nu_grav	= nugrav_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN;
	fprintf(ofv,"%e\t%e\t%e\t%e\t%e\n",t,r_hstar,nu_trans,nu_coll,nu_grav);
	fclose(ofv);
}

void problem_inloop(){
}

