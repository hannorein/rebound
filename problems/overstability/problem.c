#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
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
double buffer_zone = 1000;
char	dirname[256];
double	tau = 1.3;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 4e-3*2.*M_PI/OMEGA;
	
	sprintf(dirname,"out");	
	char* argument;
	if ((argument = input_get_argument(argc, argv, "tau"))){
		tau = atof(argument);
		sprintf(dirname,"%s_tau%e",dirname,tau);
	}
	
	if ((argument = input_get_argument(argc, argv, "eps"))){
		coefficient_of_restitution 	= atof(argument);
		sprintf(dirname,"%s_eps%e",dirname,coefficient_of_restitution);
	}else{
		coefficient_of_restitution 	= 0.5;
	}

	// reduced by 2
	root_nx 	= 500;
	root_ny 	= 1;
	root_nz 	= 8;

	double particle_r 	= 1;
	boxsize 		= 15.534876239824986;
	nghostx = 1; 	nghosty = 1; 	nghostz = 0;
	tmax = 2.*M_PI*10000.;
	init_box();
	minimum_collision_velocity = OMEGA*particle_r*0.005;



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

#ifdef MPI
	if (mpi_id==0){
#endif // MPI

	char tmp[512];
	sprintf(tmp, "rm -rf %s", dirname);
	system(tmp);
	sprintf(tmp, "mkdir %s", dirname);
	system(tmp);
#ifdef MPI
	}
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
}

int x_counter = 0;
void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (output_check(2.*M_PI)){
		char tmp[1024];
		sprintf(tmp,"%s/x.bin_%09d",dirname,x_counter);
		output_x(tmp);
		x_counter++;
	}
}

void problem_finish(){
	char tmp[256];
	sprintf(tmp,"%s/position.txt",dirname);
	output_ascii(tmp);
}
