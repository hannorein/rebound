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

extern double OMEGA;
extern double OMEGAZ;
extern double coefficient_of_restitution; 
extern double minimum_collision_velocity;
char	dirname[256];
double*	heights;
double  height;
int	heights_N = 0;
int	heights_Nmax =100;
double	tau, t_start, particle_r;
extern double minimum_collision_velocity;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 				= 1.;
	OMEGAZ 				= 3.6;
	dt				= 1e-3*2.*M_PI/OMEGA;
	tmax				= 100.*2.*M_PI;
	t_start				= tmax/2.;
	
	sprintf(dirname,"out");	
	char* argument;
	if ((argument = input_get_argument(argc, argv, "tau"))){
		tau = atof(argument);
		sprintf(dirname,"%s_tau%e",dirname,tau);
	}else{
		tau = 1.63;
	}
	
	if ((argument = input_get_argument(argc, argv, "eps"))){
		coefficient_of_restitution 	= atof(argument);
		sprintf(dirname,"%s_eps%e",dirname,coefficient_of_restitution);
	}else{
		coefficient_of_restitution 	= 0.5;
	}
	puts(dirname);

	// reduced by 2
	root_nx 	= 1;
	root_ny 	= 1;
	root_nz 	= 4;

	particle_r 	= 1;
	boxsize 	= 25;
	nghostx = 1; 	nghosty = 1; 	nghostz = 0;
	init_box();
	minimum_collision_velocity = OMEGA*particle_r*0.005;


	// Initial conditions
	double _N = tau * boxsize_x * boxsize_y/(M_PI*particle_r *particle_r);
	while (N<_N){
		struct particle p;
		p.x 	= boxsize_x*tools_uniform(-0.5,0.5);
		p.y 	= boxsize_y*tools_uniform(-0.5,0.5);
		p.z 	= particle_r*tools_uniform(-5,5);
		p.vx 	= 0;
		p.vy 	= -1.5*p.x*OMEGA;
		p.vz 	= 0;
		p.ax 	= 0; p.ay 	= 0; p.az 	= 0;
		p.m 	= 1.;
		p.r 	= particle_r;
		particles_add(p);
	}
	heights_Nmax*=N;
	heights = calloc(sizeof(double),heights_Nmax);


	char tmp[512];
	sprintf(tmp, "rm -rf %s", dirname);
	system(tmp);
	sprintf(tmp, "mkdir %s", dirname);
	system(tmp);
}



void problem_inloop(){
}

double 	veldisp_sum_xx  = 0;
double 	veldisp_sum_yy  = 0;
double 	veldisp_sum_zz  = 0;
double 	veldisp_sum_xy  = 0;
double	veldisp_sum   	= 0;
int 	veldisp_N 	= 0;
void veldisp_meassure(){
	double xx = 0;
	double yy = 0;
	double xy = 0;
	double zz = 0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		double vx = p.vx;
		double vy = p.vy+1.5*OMEGA*p.x;
		double vz = p.vz;
		xx += vx*vx;
		yy += vy*vy;
		zz += vz*vz;
		xy += vx*vy;
	}
	veldisp_sum_xx += xx/(double)N;
	veldisp_sum_yy += yy/(double)N;
	veldisp_sum_zz += zz/(double)N;
	veldisp_sum_xy += xy/(double)N;
	veldisp_sum 	= veldisp_sum_xx + veldisp_sum_yy + veldisp_sum_zz;
	veldisp_N++;

}

extern long collisions_Nlog;
extern double collisions_plog;

// This function calculates the translational part of the viscosity
double calculate_viscosity_trans_xy(){
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
double calculate_viscosity_trans_xx(){
	double nu = 0;
	double mtot = 0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		mtot	+= p.m;
		nu	+= p.m*p.vx*p.vx;
	}
	return 2./(3.*OMEGA*mtot)*nu;
}
double calculate_viscosity_trans_yy(){
	double nu = 0;
	double mtot = 0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		mtot	+= p.m;
		double vyr = p.vy + 1.5 * OMEGA * p.x;
		nu	+= p.m*vyr*vyr;
	}
	return 2./(3.*OMEGA*mtot)*nu;
}

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

double nucoll_tot=0, nutrans_yy_tot=0, nutrans_xx_tot=0, nutrans_xy_tot=0; int nuN=0;
void viscosity_meassure(){
	// Calculate viscosity
	nutrans_xy_tot 	+= calculate_viscosity_trans_xy();
	nutrans_xx_tot 	+= calculate_viscosity_trans_xx();
	nutrans_yy_tot 	+= calculate_viscosity_trans_yy();
	nucoll_tot	 = calculate_viscosity_coll();
	nuN++;
}

double fillingfactor = 0;
int fillingfactor_N = 0;

void fillingfactor_measure(){
	double area = 0;
	double H2 = 0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		double a2 = p.r*p.r - p.z*p.z;
		H2 += p.z*p.z;
		if (a2>0.){
			area += M_PI*a2;
		}
	}
	double totalarea = boxsize_x*boxsize_y;
	fillingfactor += area/totalarea;
	fillingfactor_N++;
}

void output_everything(){
	char tmp[256];
	sprintf(tmp,"%s/data.txt",dirname);
	FILE* of = fopen(tmp,"w"); 
	fprintf(of,"%e\t",tau);										//1	tau
	fprintf(of,"%e\t",coefficient_of_restitution);							//2	coeff
	fprintf(of,"%e\t",nutrans_xx_tot/(double)nuN/((double)N/(boxsize_x*boxsize_y)));		//3
	fprintf(of,"%e\t",nutrans_yy_tot/(double)nuN/((double)N/(boxsize_x*boxsize_y)));		//4
	fprintf(of,"%e\t",nutrans_xy_tot/(double)nuN/((double)N/(boxsize_x*boxsize_y)));		//5
	fprintf(of,"%e\t",veldisp_sum   /(double)veldisp_N);						//6	total veldisp^2
	fprintf(of,"%e\t",nucoll_tot/(double)nuN/((double)N/(boxsize_x*boxsize_y)));			//7
	fprintf(of,"%e\t",((double)collisions_Nlog)/(double)N/(double)(t-t_start));			//8	collision frequency
	fprintf(of,"%e\t",veldisp_sum_xx /(double)veldisp_N);						//9	vdisp_xx^2
	fprintf(of,"%e\t",veldisp_sum_yy /(double)veldisp_N);						//10	vdisp_yy^2
	fprintf(of,"%e\t",veldisp_sum_zz /(double)veldisp_N);						//11	vdisp_zz^2
	fprintf(of,"%e\t",veldisp_sum_xy /(double)veldisp_N);						//12	vdisp_xy^2
	fprintf(of,"%e\t",fillingfactor  /(double)fillingfactor_N);					//13	FF_0
	fprintf(of,"%e\t",0.);										//14	H scale height
	fprintf(of,"%e\t",fillingfactor  /(double)fillingfactor_N/(4./3.*M_PI*powf(particle_r,3.)));	//15	n volumetric number density
	fprintf(of,"%e\t",(double)N  /(boxsize_x*boxsize_y) );						//16	N surface number density
	fprintf(of,"\n");
	fclose(of);
}

void height_measure(){
	int i = 0;
	int didchange = 0;
	while(heights_N < heights_Nmax && i<N){
		heights[heights_N] = particles[i].z; 
		heights_N++;
		i++;
		didchange = 1;
	}
	if (didchange && heights_N==heights_Nmax){	
		char tmp[256];
		sprintf(tmp,"%s/heights.txt",dirname);
		FILE* file = fopen(tmp,"w");
		for (int i=0;i<heights_Nmax;i++){
			fprintf(file,"%f\n",heights[i]);
		}
		fclose(file);
	}
}


double reset_counters = 0;
void problem_output(){
	if (output_check(2.*M_PI)){
		output_timing();
	}
	if (t>t_start){
		if (reset_counters==0){
			reset_counters=1;
			collisions_Nlog = 0;
			collisions_plog = 0;
		}
		if (output_check(1e-2*2.*M_PI)){
			veldisp_meassure();
			viscosity_meassure();
			fillingfactor_measure();
			height_measure();
		}
	}
}

void problem_finish(){
	output_everything();
	char tmp[256];
	sprintf(tmp,"%s/position.txt",dirname);
	output_ascii(tmp);
}
