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
#include <string.h>
#ifdef OPENMP
#include <omp.h>
#endif
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "raytracing.h"
#include "collisions.h"
#include "collision_resolve.h"
#include "tools.h"
#include "input.h"

extern double OMEGA;

void input_append_input_arguments_with_int(const char* argument, int value);
void problem_find_restartfile();


extern double opening_angle2;
extern int output_logfile_first;

double coefficient_of_restitution_bridges(double v){ // assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void problem_init(int argc, char* argv[]){
	// Setup constants
	opening_angle2	= .5;
	G 				= 6.67428e-11;		// N m^2 / kg^2 
	double a			= 391e3;		// m 
	double chariklo_m		= 1e19; 		// kg
	OMEGA 				= sqrt(G*chariklo_m/(a*a*a));	// 1/s
	softening 			= 0.1;			// m
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
	tmax				= 10.*2.*M_PI/OMEGA+dt;

	// Setup domain dimensions
	int _root_n = 1;

	root_nx = _root_n;
	root_ny = _root_n;
	root_nz = 1;
	nghostx = 2;
	nghosty = 2;
	nghostz = 0;
	boxsize = input_get_double(argc,argv,"boxsize",25)/(double)root_nx;
	init_box();
	
	// Setup particle and disk properties
	double sigma 			= input_get_double(argc,argv,"sigma",400); 			// kg/m^2
	double rho			= input_get_double(argc,argv,"rho",250);			// kg/m^3
	double particle_radius_min 	= input_get_double(argc,argv,"rmin",.3);				// m
	double particle_radius_max 	= input_get_double(argc,argv,"rmax",3);				// m
	double particle_radius_slope 	= -3;	
	minimum_collision_velocity 	= 0.05*particle_radius_min*OMEGA;
	coefficient_of_restitution_for_velocity	= coefficient_of_restitution_bridges;

	
	struct 	aabb bb	= {	.xmin = -boxsize_x/2., 
				.xmax =  boxsize_x/2., 
				.ymin = -boxsize_y/2., 
				.ymax =  boxsize_y/2., 
				.zmin = -boxsize_z/2., 
				.zmax =  boxsize_z/2.	};
	long	_N;
	if (particle_radius_min==particle_radius_max){
		_N = round(sigma*boxsize_x*boxsize_y/(4./3.*M_PI*rho*  pow(particle_radius_max,3.)));
	}else{
		_N = round(sigma*boxsize_x*boxsize_y/(4./3.*M_PI*rho* (pow(particle_radius_max,4.+particle_radius_slope) - pow(particle_radius_min,4.+particle_radius_slope)) / (pow(particle_radius_max,1.+particle_radius_slope) - pow(particle_radius_min,1.+particle_radius_slope)) * (1.+particle_radius_slope)/(4.+particle_radius_slope)));
	}

#ifdef OPENMP
	input_append_input_arguments_with_int("openmp",omp_get_max_threads());
#endif // OPENMP

	if (output_check_directory()==1){
		// Restart simulation
		if( access( "binary.txt", F_OK ) == -1 ) {
			printf("Cannot find a restart file in directory.\n");
			exit(-1);
		}
		input_binary("binary.txt");	
		output_logfile_first = 0;
		output_double("Restarted at time",t);
	}else{
		// Start from scratch and prepare directory
		output_prepare_directory();
		double mass = 0;
		double area = 0;
		for(int i=0;i<_N;i++){
			struct particle pt;
			do{
				pt.z 		= particle_radius_max*tools_normal(2.);	
			}while(fabs(pt.z)>=boxsize_z/2.);
			pt.x 		= tools_uniform(bb.xmin,bb.xmax);
			pt.y 		= tools_uniform(bb.ymin,bb.ymax);
			pt.vx 		= 0;
			pt.vy 		= -1.5*pt.x*OMEGA;
			pt.vz 		= 0;
			pt.ax 		= 0;
			pt.ay 		= 0;
			pt.az 		= 0;
			pt.r 		= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
			pt.m 		= rho*4./3.*M_PI*pt.r*pt.r*pt.r;
			mass		+= pt.m;
			area		+= M_PI*pt.r*pt.r;

			
			particles_add(pt);
		}
		output_double("boxsize [m]",boxsize);
		output_int("root_nx",root_nx);
		output_int("root_ny",root_ny);
		output_int("root_nz",root_nz);
		output_double("tau_geom",area/(boxsize*boxsize));
		output_double("sigma [kg/m^2]",sigma);
		output_double("rho [kg/m^3]",rho);
		output_double("particle_radius_min [m]",particle_radius_min);
		output_double("particle_radius_max [m]",particle_radius_max);
		output_double("particle_radius_slope",particle_radius_slope);
		output_double("OMEGA [1/s]",OMEGA);
		output_double("lambda_crit [m]",4.*M_PI*M_PI*G*sigma/OMEGA/OMEGA);
		output_double("length [m]",boxsize*(double)root_nx);
		output_double("length/lambda_crit",boxsize*(double)root_nx/(4.*M_PI*M_PI*G*sigma/OMEGA/OMEGA));
		output_int("N",_N);
		output_double("tmax [orbits]",tmax/(2.*M_PI/OMEGA));
		output_int("number of timesteps",ceil(tmax/dt));
		system("cat config.log");
	}
}


void problem_inloop(){
}

void problem_output(){
	if (output_check(10.*dt)){
		int N_rays 	= 1000; 
		double B 	= M_PI_2; 	// elevation angle of observer (here: normal to midplane)
		double sun_B 	= M_PI_2; 	// elevation angle of light source (here: normal to midplane)
		double phi	= 0;		// azimuthal angle of observer (here: doesn't matter as B=M_PI_2)
		double sun_phi	= 0;		// azimuthal angle of light source (here: doesn't matter as B_sun=M_PI_2)
		double flux;			// this will be set to the flux (0..1)
		double opacity;			// this will be set to the opacity (0..1)
		tree_raytrace(N_rays, B, phi, sun_B, sun_phi, &flux, &opacity);

		printf("\n flux: %0.6f \t opacity: %0.6f\n",flux,opacity);
	}
	if (output_check(10.*dt)){
		output_timing();
	}
	if (output_check(.2*M_PI/OMEGA)){
		output_binary("binary.txt");
	}
}

void problem_finish(){
}
