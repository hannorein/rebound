#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <fftw3.h>
#ifdef OPENMP
#include <omp.h>
#endif //OPENMP
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
extern double collisions_plog;
extern long   collisions_Nlog;
extern double (*coefficient_of_restitution_for_velocity)(double); 
void input_append_input_arguments_with_double(const char* argument, double value);
double coefficient_of_restitution_bridges(double v); 
double buffer_zone = 0;
double	tau;
double particle_r;
double particle_min_r;
double particle_max_r;
const int fft_N = 1024;


void problem_init(int argc, char* argv[]){
	int id = input_get_int(argc,argv,"id",0)-1;
	if (id<0){
		tau 					= input_get_double(argc, argv, "tau",1.64);
		coefficient_of_restitution 		= input_get_double(argc, argv, "eps",0.5);
	}else{
		// Setup a grid of 16x16 simulations
		int ntau = 16;
		int neps = 16;
		int itau = id%neps;
		int ieps = id/neps;
		tau 					= 0.1+2.5*(double)itau/(double)(ntau-1);
		coefficient_of_restitution 		= 0.0+0.9*(double)ieps/(double)(neps-1);
		input_append_input_arguments_with_double("tau",tau);
		input_append_input_arguments_with_double("eps",coefficient_of_restitution);
	}
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "dtorb",4e-3);
	tmax 				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "tmaxorb",10000);
	
	buffer_zone			= input_get_double(argc, argv, "buffer_zone",0);
	double delta_tau		= input_get_double(argc, argv, "delta_tau",0);
	if (input_get_int(argc,argv,"bridges",0)){
		coefficient_of_restitution_for_velocity = coefficient_of_restitution_bridges;
	}
	
	root_nx 			= input_get_int(argc,argv,"root_nx",64);
	root_ny 			= 1;
	root_nz 			= 16;
	nghostx 			= 1; 	
	nghosty 			= 1; 	
	nghostz 			= 0;

	particle_min_r 			= input_get_double(argc, argv, "particle_min_r",1);
	particle_max_r			= input_get_double(argc, argv, "particle_max_r",1);
	particle_r 			= (particle_max_r+particle_min_r)/2.;
	boxsize 			= input_get_double(argc, argv, "boxsize",15.534876239824986);
	init_box();
	system("mkdir out");
	chdir("out");
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

	double area_total = tau * (xmax-xmin) * boxsize_y;
	double area_added = 0;

	while (area_added<area_total&&xmax>xmin){
		struct particle p;
		double x,prob;
		do{
			x 	= tools_uniform(xmin,xmax);
			prob 	= tools_uniform(0,tau+delta_tau);
		}while(prob> x/(xmax-xmin)*2.*delta_tau+tau);
		p.x 	= x;
		p.y 	= ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z 	= 10.0*((double)rand()/(double)RAND_MAX-0.5)*particle_r;
		p.vx 	= 0;
		p.vy 	= -1.5*p.x*OMEGA;
		p.vz 	= 0;
		p.ax 	= 0; p.ay 	= 0; p.az 	= 0;
		p.r 	= tools_powerlaw(particle_min_r, particle_max_r,-3);
		p.m 	= p.r*p.r*p.r;
		area_added += p.r*p.r*M_PI;
		particles_add(p);
	}

#ifdef MPI
	if (mpi_id==0){
	output_int("mpi_num",mpi_num);
#endif // MPI
	output_int("N",N);
	output_double("dt",dt);
	output_double("tau",tau);
	output_double("delta_tau",delta_tau);
	output_double("coefficient_of_restitution",coefficient_of_restitution);
	output_double("boxsize",boxsize);
	output_double("boxsize_x",boxsize_x);
	output_double("boxsize_y",boxsize_y);
	output_double("boxsize_z",boxsize_z);
	output_double("particle_min_r",particle_min_r);
	output_double("particle_max_r",particle_max_r);
	output_double("buffer_zone",buffer_zone);
	output_double("fft_lamda_min",boxsize_x/(double)fft_N);

#ifdef MPI
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif // MPI
}

double coefficient_of_restitution_bridges(double v){
	v *= 0.00013143527;	
	// assumes v in units of [m/s]
	double eps = 0.32*pow(fabs(v)*100.,-0.234);
	if (eps>1) eps=1;
	if (eps<0) eps=0;
	return eps;
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
		double x = tools_uniform(-boxsize_x/2.+particle_max_r,boxsize_x/2.-particle_max_r); // random position
		double y = tools_uniform(-boxsize_y/2.+particle_max_r,boxsize_y/2.-particle_max_r); // random position
		int i1 = particle_close_to(x-particle_max_r)-1;
		i1 = i1<0?0:i1;
		int i2 = particle_close_to(x+particle_max_r)+1;
		i2 = i2>=N?N-1:i2;
		int tau = 0;
		for(int j=i1; j<i2; j++){
			struct particle p = particles[j];
			double distance2 = (p.y-y)*(p.y-y)+(p.x-x)*(p.x-x);
			if (distance2<p.r*p.r){
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

void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (output_check(12.242342345*2.*M_PI)){
		output_append_tau("tau.txt");
	}
}

void problem_finish(){
	output_ascii("position_end.txt");
}
