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
double buffer_zone = 0;
double	tau;
double particle_r;
const int fft_N = 1024;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "dtorb",4e-3);
	tmax 				= 2.*M_PI/OMEGA*input_get_double(argc, argv, "tmaxorb",500.);
	
	buffer_zone			= input_get_double(argc, argv, "buffer_zone",0);
	tau 				= input_get_double(argc, argv, "tau",1.64);
	double delta_tau		= input_get_double(argc, argv, "delta_tau",0);
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
		double x,prob;
		do{
			x 	= tools_uniform(xmin,xmax);
			prob 	= tools_uniform(0,tau+delta_tau/2.);
		}while(prob> x/(xmax-xmin)*2.*delta_tau+tau);
		p.x 	= x;
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
	output_double("particle_r",particle_r);
	output_double("buffer_zone",buffer_zone);
	output_double("fft_lamda_min",boxsize_x/(double)fft_N);

#ifdef MPI
	}
	MPI_Barrier(MPI_COMM_WORLD);
#endif // MPI
}

void output_append_energy(char* filename){
	double energy = 0; 
#pragma omp parallel for reduction(+:energy)
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		const double vy = p.vy + 1.5*p.x*OMEGA;
		energy += 0.5*(vy*vy + p.vx*p.vx);
	}
	FILE* of = fopen(filename,"a+"); 
	fprintf(of,"%e\t%e\n",t,energy/(double)N);
	fclose(of);
}

void output_x(char* filename1){
	char filename[1024];
	sprintf(filename,"%s_%015.3f",filename1,t);
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

fftw_complex* fft_in 	= NULL;
double* fft_in_N 	= NULL;
fftw_complex* fft_out	= NULL;
fftw_plan fft_p		= NULL;

void fft_empty_data(){
#pragma omp parallel for 
	for (int i=0;i<fft_N;i++){
		fft_in[i][0] 	= 0;
		fft_in[i][1] 	= 0;
		fft_in_N[i] 	= 0;
	}
}
struct fft_output{
	double power;
	double lambda;	
};

int compare_fft_output (const void * a, const void * b){
	const double diff = ((struct fft_output*)a)->power - ((struct fft_output*)b)->power;
	if (diff < 0) return 1;
	if (diff > 0) return -1;
	return 0;
}

void fft_write_to_file(char* filename){
	struct fft_output* fft_outputs = malloc(sizeof(struct fft_output)*(fft_N/2-1));
	for (int i=1;i<fft_N/2;i++){
		fft_outputs[i-1].power	= fft_out[i][0]*fft_out[i][0] + fft_out[i][1]*fft_out[i][1];
		fft_outputs[i-1].lambda	= boxsize_x/(double)i;
	}
	qsort (fft_outputs, fft_N/2-1, sizeof(struct fft_output), compare_fft_output);
	FILE* of = fopen(filename,"a+"); 
	// Output ten highest power modes
	for (int i=0;i<10;i++){
		fprintf(of,"%e\t%e\t%e\n",t,fft_outputs[i].lambda,fft_outputs[i].power);
	}
	fclose(of);
}

void output_fft(){
	if (fft_in==NULL){
		fft_in_N 		= calloc(fft_N,sizeof(double));
		fft_in 			= (fftw_complex*)malloc(sizeof(double)*2*fft_N);
		fft_out 		= (fftw_complex*)malloc(sizeof(double)*2*fft_N);
#if OPENMP
		fftw_plan_with_nthreads(omp_get_num_threads());
#endif // OPENMP
		fft_p			= fftw_plan_dft_1d(fft_N, fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
	}
	const double sqrt_fft_N = sqrt(fft_N);

	// ******* Density ********
	fft_empty_data();
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		int j = ((int)(floor((p.x + boxsize_x/2.)/boxsize_x*(double)fft_N))+fft_N)%fft_N;
		fft_in[j][0] += 1;
	}
#pragma omp parallel for 
	for (int i=0;i<fft_N;i++){
		fft_in[i][0] *= 1./(double)N * (double)(fft_N)/sqrt_fft_N;
	}
	fftw_execute(fft_p);
	fft_write_to_file("fft_density.txt");

	// ******* Velocity ********
	/*
	fft_empty_data();
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		int j = ((int)(floor((p.x + boxsize_x/2.)/boxsize_x*(double)fft_N))+fft_N)%fft_N;
		fft_in[j][0] += sqrt(p.vx*p.vx + (p.vy+1.5*OMEGA*p.x)*(p.vy+1.5*OMEGA*p.x) + p.vz*p.vz);
		fft_in_N[j]  += 1.0;
	}
#pragma omp parallel for 
	for (int i=0;i<fft_N;i++){
		if (fft_in_N[i]>0){
			fft_in[i][0] /= fft_in_N[i] * sqrt_fft_N;
		}
	}
	fftw_execute(fft_p);
	fft_write_to_file("fft_velocity.txt");
	*/
}

extern void collisions_sweep_shellsort_particles();
void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (t>2.*M_PI*400&&t<2.*M_PI*500.&&output_check(.02*M_PI)){
		output_x("xfine.bin");
	}
	if (output_check(20.*M_PI)){
		output_x("x.bin");
	}
	if (output_check(2.*M_PI)){
		output_append_velocity_dispersion("vdisp.txt");
		output_append_energy("energy.txt");
	}
	if (output_check(50.*2.*M_PI)){
		char filename[1024];
		sprintf(filename,"ascii_%05d.txt",(int)round(t/2./M_PI));
		output_ascii(filename);
	}
	if (output_check(12.242342345*2.*M_PI)){
		output_append_tau("tau.txt");
		output_fft();
	}
}

void problem_finish(){
	output_ascii("position_end.txt");
}
