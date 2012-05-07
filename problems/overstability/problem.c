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
double particle_r;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "dtorb",4e-3);
	tmax 				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "tmaxorb",10000.);
	
	tau 				= input_get_double(argc, argv, "tau",1.64);
	coefficient_of_restitution 	= input_get_double(argc, argv, "eps",0.5);
	root_nx 			= input_get_int(argc,argv,"root_nx",500);
	root_ny 			= 1;
	root_nz 			= 8;
	nghostx 			= 1; 	
	nghosty 			= 1; 	
	nghostz 			= 0;

	particle_r 			= input_get_double(argc, argv, "particle_r",1);
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
int particle_close_to(double x){
	int i_min = 0;
	int i_max = N-1;
	while(i_min < i_max-1){
		int i_middle = (i_min+i_max)/2;
		if (x <= particles[i_middle].x){
			i_max = i_middle;
		}else{
			i_min = i_middle;
		}
	}
	return i_min;
}

void output_append_tau(char* filename){
	FILE* of = fopen(filename,"a"); 
	if (of==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}
	long _N = 100000; // Number of Monte Carlo samples.
	double tau_geometric 	= 0;
	double tau_photometric 	= 0;
#pragma omp parallel for reduction(+:tau_geometric) reduction(+:tau_photometric)
	for (int i=0;i<_N;i++){
		double x = tools_uniform(-boxsize_x/2.+particle_r,boxsize_x/2.-particle_r); // random position
		double y = tools_uniform(-boxsize_y/2.+particle_r,boxsize_y/2.-particle_r); // random position
		int i1 = particle_close_to(x-particle_r)-1;
		i1 = i1<0?0:i1;
		int i2 = particle_close_to(x+particle_r)+1;
		i2 = i2>=N?N-1:i2;
		int tau = 0;
		for(int j=i1; j<i2; j++){
			struct particle p = particles[j];
			double distance2 = (p.y-y)*(p.y-y)+(p.x-x)*(p.x-x);
			if (distance2<particle_r*particle_r){
				tau+=1;
			}
		}
		tau_photometric	+= tau?1.0:0.0;
		tau_geometric	+= (double)tau;
	}
	tau_geometric	/= (double)_N;
	tau_photometric	/= (double)_N;
	fprintf(of,"%e\t%e\t%e\n",t,tau_geometric, tau_photometric);
	fclose(of);

}

extern void collisions_sweep_insertionsort_particles();
void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (output_check(2.*M_PI)){
		char tmp[1024];
		sprintf(tmp,"x.bin_%015.3f",t);
		output_x(tmp);
		output_append_velocity_dispersion("vdisp.txt");
	}
	if (output_check(.2*M_PI)){
		// Optical depth
		collisions_sweep_insertionsort_particles();
		output_append_tau("tau.txt");
	}
}

void problem_finish(){
	output_ascii("position_end.txt");
}
