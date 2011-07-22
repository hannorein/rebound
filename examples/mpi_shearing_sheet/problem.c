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

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 


void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 	= 0.00013143527;		// 1/s
	G 	= 6.67428e-11;			// m^3 / kg s^2
	softening = 0.1;			// m
	dt = 1e-3*2.*M_PI/OMEGA;		// s
	nghostx = 3; nghosty = 3; nghostz = 0; 	// Use three ghost rings
	double surfacedensity 	= 400; 		// kg/m^2
	double particle_density	= 400;		// g/cm^3
	double particle_radius 	= 5;		// m
	double particle_mass 	= particle_density*4./3.*M_PI*pow(particle_radius,3); 	// kg
	if (argc>1){				// Try to read boxsize from command line
		boxsize = atof(argv[1]);
	}else{
		boxsize = 100;
	}
	root_nx = 2; root_ny = 2; root_nz = 2;
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius*OMEGA*0.01;  // small fraction of the shear
	// Setup particle structures
	int _N = round(surfacedensity*boxsize*boxsize/particle_mass);
	init_box();
	// Initial conditions
#ifdef MPI
	// Only initialise particles on master. This should also be parallised, but obviously depends on the problem.
	if (mpi_id==0){
#endif
		for (int i =0;i<_N;i++){
			struct particle pt;
			double vrand = 0.01*OMEGA*((double)rand()/(double)RAND_MAX-0.5);
			double phirand = 2.*M_PI*((double)rand()/(double)RAND_MAX-0.5);
			pt.x 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
			pt.y 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
			pt.z 		= 10.*((double)rand()/(double)RAND_MAX-0.5);
			pt.vx 		= 0;
			pt.vy 		= -1.5*pt.x*OMEGA+2.*vrand*cos(phirand);
			pt.vz 		= vrand*sin(phirand);
			pt.ax 		= 0;
			pt.ay 		= 0;
			pt.az 		= 0;
			pt.m 		= particle_mass;
#ifndef COLLISIONS_NONE
			pt.r 		= particle_radius;
#endif
			particles_add(pt);
		}
#ifdef MPI
	}
#endif
	printf("%d particles initialized on root node\n",_N);
}

double coefficient_of_restitution_bridges(double v){
	// v in [m/s]
	return 0.34*pow(fabs(v)*100.,-0.234);
}

void problem_inloop(){

}

void problem_output(){
	if (output_check(2.*M_PI/OMEGA)){
		output_timing();
	}
}

void problem_finish(){
}
