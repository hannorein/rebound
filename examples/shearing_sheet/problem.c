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

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 

extern double opening_angle2;

void problem_init(int argc, char* argv[]){
	// Setup constants
#ifdef GRAVITY_TREE
	opening_angle2	= .5;
#endif
	OMEGA 		= 0.00013143527;		// 1/s
	G 		= 6.67428e-11;			// N / (1e-5 kg)^2 m^2
	softening 	= 0.1;				// m
	dt 		= 1e-3*2.*M_PI/OMEGA;		// s
	root_nx = 2; root_ny = 2; root_nz = 1;
	nghostx = 2; nghosty = 2; nghostz = 0; 		// Use three ghost rings
	double surfacedensity 	= 400; 			// kg/m^2
	double particle_density	= 400;			// kg/m^3
	double particle_radius_min = 1;			// m
	double particle_radius_max = 4;			// m
	double particle_radius_slope = -3;	
	if (argc>1){					// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}else{
		boxsize = 100;
	}
	printf("Toomre wavelength: %f\n",2.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*G);
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear
	// Setup particle structures
	init_box();
	double total_mass = surfacedensity*boxsize_x*boxsize_y;
	// Initial conditions
#ifdef MPI
	// Only initialise particles on master. This should also be parallised, but obviously depends on the problem.
	if (mpi_id==0){
#endif
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
#ifdef MPI
	printf("%d particles initialized on root node.\n",_N);
	}
#endif
}

double coefficient_of_restitution_bridges(double v){
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void problem_inloop(){

}

void problem_output(){
#ifdef LIBPNG
	if (output_check(1e-3*2.*M_PI/OMEGA)){
		output_png("png/");
	}
#endif //LIBPNG
	/*
	if (output_check(1e-1*2.*M_PI/OMEGA)){
		output_timing();
		output_append_velocity_dispersion("veldisp.txt");
	}
	if (output_check(2.*M_PI/OMEGA)){
		output_ascii("position.txt");
	}
	*/
}

void problem_finish(){
}
