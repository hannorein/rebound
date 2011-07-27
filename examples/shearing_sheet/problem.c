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
		boxsize = 200;
	}
	// Use Bridges et al coefficient of restitution.
	coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	minimum_collision_velocity = particle_radius*OMEGA*0.1;  // small fraction of the shear
	// Setup particle structures
	init_box();
	int _N = (int)round(surfacedensity*boxsize*boxsize/particle_mass);
	// Initial conditions
	for (int i =0;i<_N;i++){
		struct particle p;
		double vrand = 0.0*OMEGA*((double)rand()/(double)RAND_MAX-0.5);
		p.x 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y 		= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z 		= 0;
		p.vx 		= 0;
		p.vy 		= -1.5*p.x*OMEGA;
		p.vz 		= 0;;
		p.ax 		= 0;
		p.ay 		= 0;
		p.az 		= 0;
		p.m 		= particle_mass;
		p.r 		= particle_radius;
		particles_add(p);
	}
}

double coefficient_of_restitution_bridges(double v){
	// v in [m/s]
	double eps = 0.34*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
}

void problem_inloop(){

}

void problem_output(){
	if (output_check(1e-1*2.*M_PI/OMEGA)){
		output_append_velocity_dispersion("veldisp.txt");
	}
}

void problem_finish(){
}
