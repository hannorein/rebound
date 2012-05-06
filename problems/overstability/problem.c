#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "main.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"
#include "input.h"
#include "output.h"
#include "collisions.h"
#include "communication_mpi.h"

extern double OMEGA;
extern double OMEGAZ;
extern double coefficient_of_restitution; 
extern double minimum_collision_velocity;
double buffer_zone = 0;
double	tau;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "dtorb",4e-3);
	tmax 				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "tmaxorb",10000.);
	
	tau 				= input_get_double(argc, argv, "tau",1.3);
	coefficient_of_restitution 	= input_get_double(argc, argv, "eps",0.5);
	root_nx 			= input_get_int(argc,argv,"root_nx",500);
	root_ny 			= 1;
	root_nz 			= 8;
	nghostx 			= 1; 	
	nghosty 			= 1; 	
	nghostz 			= 0;

	double particle_r 		= input_get_double(argc, argv, "particle_r",1);
	boxsize 			= input_get_double(argc, argv, "boxsize",15.534876239824986);
	init_box();
	output_prepare_directory();
	
	// Initial conditions
#ifdef MPI
	double xmin = -boxsize_x/2.+(double)mpi_id*boxsize_x/(double)mpi_num;
	double xmax = -boxsize_x/2.+(double)(mpi_id+1)*boxsize_x/(double)mpi_num;
#else // MPI
	double xmin = -boxsize_x/2.;
	double xmax = boxsize_x/2.;
#endif // MPI
	if (xmin<-boxsize_x/2.+buffer_zone){
		xmin = -boxsize_x/2.+buffer_zone;
	}
	if (xmax<-boxsize_x/2.+buffer_zone){
		xmax = -boxsize_x/2.+buffer_zone;
	}
	if (xmin>boxsize_x/2.-buffer_zone){
		xmin = boxsize_x/2.-buffer_zone;
	}
	if (xmax>boxsize_x/2.-buffer_zone){
		xmax = boxsize_x/2.-buffer_zone;
	}

	double _N = tau * (xmax-xmin) * boxsize_y/(M_PI*particle_r *particle_r);
	while (N<_N&&xmax>xmin){
		struct particle p;
		p.x 	= tools_uniform(xmin,xmax);
		p.y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z 	= 10.0*((double)rand()/(double)RAND_MAX-0.5)*particle_r;
		p.vx 	= 0;
		p.vy 	= -1.5*p.x*OMEGA;
		p.vz 	= 0;
		p.ax 	= 0; p.ay 	= 0; p.az 	= 0;
		p.m 	= 1.;
		p.r 	= particle_r;
		particles_add(p);
	}

	output_int("N",N);
	output_double("dt",dt);
	output_double("tau",tau);
	output_double("coefficient_of_restitution",coefficient_of_restitution);
	output_double("boxsize",boxsize);
	output_double("boxsize_x",boxsize_x);
	output_double("boxsize_y",boxsize_y);
	output_double("boxsize_z",boxsize_z);
	output_double("particle_r",particle_r);

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif // MPI
}


void output_x(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"wb"); 
#else // MPI
	FILE* of = fopen(filename,"wb"); 
#endif // MPI
	if (of==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fwrite(&(p.x),sizeof(double),1,of);
	}
	fclose(of);
}

void problem_inloop(){
	//Slowly turn on the minimum collision velocity.
	const double t_init = 1./OMEGA;
	minimum_collision_velocity 	= 0.5/dt;
	if (t<t_init){
		minimum_collision_velocity 	*= t/t_init;
	}
}

void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (output_check(2.*M_PI)){
		char tmp[1024];
		sprintf(tmp,"x.bin_%09.3f",t);
		output_x(tmp);
		output_append_velocity_dispersion("vdisp.txt");
	}
}

void problem_finish(){
	output_ascii("position_end.txt");
}
