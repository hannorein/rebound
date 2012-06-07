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
#include "main.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "communication_mpi.h"
#include "tree.h"
#include "collisions.h"
#include "tools.h"
#include "input.h"

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;

extern double (*coefficient_of_restitution_for_velocity)(double); 
double coefficient_of_restitution_bridges(double v); 
void input_append_input_arguments_with_int(const char* argument, int value);
void problem_find_restartfile();
int position_id=0;	// output number


extern double opening_angle2;
extern int output_logfile_first;

void problem_init(int argc, char* argv[]){
	// Setup constants
#ifdef GRAVITY_TREE
	opening_angle2	= .5;
#endif
	OMEGA 				= 0.0001948170;		// 1/s corresponds to a=100e6 m
	G 				= 6.67428e-11;		// N m^2 / kg^2 
	softening 			= 0.1;			// m
	dt 				= 1e-3*2.*M_PI/OMEGA;	// s
	tmax				= 10.*2.*M_PI/OMEGA;

	// Setup domain dimensions
	int _root_n = 1;
#ifdef MPI
	do{
		_root_n++;
	}while(_root_n*_root_n <= mpi_num);
	_root_n--;
	if (_root_n*_root_n!=mpi_num){
		printf("\n\nError. mpi_num must be square of integer.\n");
		exit(-1);
	}
#endif // MPI

	root_nx = _root_n;
	root_ny = _root_n;
	root_nz = 1;
	nghostx = 2;
	nghosty = 2;
	nghostz = 0;
	boxsize = input_get_double(argc,argv,"length",500)/(double)root_nx;
	init_box();
	
	// Setup particle and disk properties
	double sigma 			= input_get_double(argc,argv,"sigma",1200); 			// kg/m^2
	double rho			= input_get_double(argc,argv,"rho",450);			// kg/m^3
	double particle_radius_min 	= input_get_double(argc,argv,"rmin",1);				// m
	double particle_radius_max 	= input_get_double(argc,argv,"rmax",1);				// m
	double particle_radius_slope 	= -3;	
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


#ifdef MPI
	input_append_input_arguments_with_int("mpinum",mpi_num);
#endif // MPI
	if (output_check_directory()==1){
		// Restart simulation
		problem_find_restartfile();
		output_logfile_first = 0;
#ifdef MPI
		if (mpi_id==0){
#endif // MPI
		output_double("Restarted at time",t);
#ifdef MPI
		}
#endif // MPI
	}else{
		// Start from scratch and prepare directory
		output_prepare_directory();
#ifdef MPI
		bb = communication_boundingbox_for_proc(mpi_id);
		_N   /= mpi_num;
		if (mpi_id==0){
#endif // MPI
			output_double("boxsize [m]",boxsize);
			output_int("root_nx",root_nx);
			output_int("root_ny",root_ny);
			output_int("root_nz",root_nz);
			output_double("tau (r_max)",_N*M_PI*particle_radius_max*particle_radius_max/(boxsize*boxsize));
			output_double("sigma [kg/m^2]",sigma);
			output_double("rho [kg/m^3]",rho);
			output_double("OMEGA [1/s]",OMEGA);
			output_double("lambda_crit [m]",4.*M_PI*M_PI*G*sigma/OMEGA/OMEGA);
			output_double("length [m]",boxsize*(double)root_nx);
			output_double("length/lambda_crit",boxsize*(double)root_nx/(4.*M_PI*M_PI*G*sigma/OMEGA/OMEGA));
			output_int("N",_N);
#ifdef MPI
			output_int("N_total",_N*mpi_num);
			output_int("mpi_num",mpi_num);
#endif // MPI
			output_double("tmax [orbits]",tmax/(2.*M_PI/OMEGA));
			output_int("number of timesteps",ceil(tmax/dt));
			system("cat config.log");
#ifdef MPI
		}
#endif // MPI
		for(int i=0;i<_N;i++){
			struct particle pt;
			do{
				pt.z 		= particle_radius_max*tools_normal(2.);	
			}while(fabs(pt.z)>=boxsize_z/2.);
#ifdef MPI
			int proc_id;
			do{
#endif // MPI
				pt.x 		= tools_uniform(bb.xmin,bb.xmax);
				pt.y 		= tools_uniform(bb.ymin,bb.ymax);
#ifdef MPI
				int rootbox = particles_get_rootbox_for_particle(pt);
				int root_n_per_node = root_n/mpi_num;
				proc_id = rootbox/root_n_per_node;
			}while(proc_id != mpi_id );
#endif // MPI
			
			pt.vx 		= 0;
			pt.vy 		= -1.5*pt.x*OMEGA;
			pt.vz 		= 0;
			pt.ax 		= 0;
			pt.ay 		= 0;
			pt.az 		= 0;
			pt.r 		= tools_powerlaw(particle_radius_min,particle_radius_max,particle_radius_slope);
			pt.m 		= rho*4./3.*M_PI*pt.r*pt.r*pt.r;
			
			particles_add(pt);
		}
	}
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif // MPI
}

void problem_find_restartfile(){
	int latest = -1;
	for (int i=0;i<1000;i++){
		char filename[256];
#ifdef MPI
		sprintf(filename,"binary_%08d.txt_0",i);
#else // MPI
		sprintf(filename,"binary_%08d.txt",i);
#endif // MPI
		if( access( filename, F_OK ) != -1 ) {
			latest = i;
		}
	}
	if (latest ==-1){
		printf("Cannot find a restart file in directory.\n");
		exit(-1);
	}
	char filename[256];
	sprintf(filename,"binary_%08d.txt",latest);
	input_binary(filename);	
}

double coefficient_of_restitution_bridges(double v){
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

void output_ascii_mod(char* filename){
#ifdef MPI
	char filename_mpi[1024];
	sprintf(filename_mpi,"%s_%d",filename,mpi_id);
	FILE* of = fopen(filename_mpi,"w"); 
#else // MPI
	FILE* of = fopen(filename,"w"); 
#endif // MPI
	if (of==NULL){
		printf("\n\nError while opening file '%s'.\n",filename);
		return;
	}
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fprintf(of,"%e\t%e\t%e\t%e\n",p.x,p.y,p.z,p.r);
	}
	fclose(of);
}
	
double a[3];
double b[3];

int tree_get_intersection(struct cell* c){
	double xmin = c->x-c->w/2.-collisions_max_r;
	double xmax = c->x+c->w/2.+collisions_max_r;
	double ymin = c->y-c->w/2.-collisions_max_r;
	double ymax = c->y+c->w/2.+collisions_max_r;
	double zmin = c->z-c->w/2.-collisions_max_r;
	double zmax = c->z+c->w/2.+collisions_max_r;

	if (b[0]!=0){
		double d,y,z;
		d = (xmin-a[0])/b[0];
		y = a[1]+d*b[1];
		z = a[2]+d*b[2];
		if (y<ymax && y>ymin && z<zmax && z>zmin) return 1;
		d = (xmax-a[0])/b[0];
		y = a[1]+d*b[1];
		z = a[2]+d*b[2];
		if (y<ymax && y>ymin && z<zmax && z>zmin) return 1;
	}
	
	if (b[1]!=0){
		double d,x,z;
		d = (ymin-a[1])/b[1];
		x = a[0]+d*b[0];
		z = a[2]+d*b[2];
		if (x<xmax && x>xmin && z<zmax && z>zmin) return 1;
		d = (ymax-a[1])/b[1];
		x = a[0]+d*b[0];
		z = a[2]+d*b[2];
		if (x<xmax && x>xmin && z<zmax && z>zmin) return 1;
	}
	
	if (b[2]!=0){
		double d,x,y;
		d = (zmin-a[2])/b[2];
		x = a[0]+d*b[0];
		y = a[1]+d*b[1];
		if (x<xmax && x>xmin && y<ymax && y>ymin) return 1;
		d = (zmax-a[2])/b[2];
		x = a[0]+d*b[0];
		y = a[1]+d*b[1];
		if (x<xmax && x>xmin && y<ymax && y>ymin) return 1;
	}

	return 0;
}

int tree_get_transparency(struct cell* c){
	if (c->pt<0){
		for(int i=0;i<8;i++){
			struct cell* o = c->oct[i];
			if (o!=NULL){
				if(tree_get_intersection(o)){
					if (tree_get_transparency(o)==1) return 1;
				}
			}
		}
		return 0;
	}else{
		struct particle p = particles[c->pt];
		double t1 = b[1]*(a[2]-p.z)-b[2]*(a[1]-p.y);
		double t2 = b[2]*(a[0]-p.x)-b[0]*(a[2]-p.z);
		double t3 = b[0]*(a[1]-p.y)-b[1]*(a[0]-p.x);
		double bot2 = b[0]*b[0]+b[1]*b[1]+b[2]*b[2];
		double distance2 = (t1*t1+t2*t2+t3*t3)/bot2;
		if (distance2<p.r*p.r) return 1;
		return 0;
	}
}

void output_transparency(double B, double phi){
	if (fabs(B-M_PI/2.)>M_PI/4.){
		printf("ERROR: Cannot calculate ray. B has to be in the range [M_PI/4.:3.*M_PI/4.].\n");
		return;
	}
	
	if (B!=M_PI/2.){
		b[0] = cos(phi)*cos(B)/sin(B);
		b[1] = sin(phi)*cos(B)/sin(B);
	}else{
		b[0] = 0;
		b[1] = 0;
	}
	a[2] = 0;
	b[2] = 1;


	int num_rays 		= 10000;

	double* _a = malloc(sizeof(double)*num_rays*3);
	int*   _t = calloc(num_rays,sizeof(int));
#ifdef MPI
	if (mpi_id==0){
#endif //MPI
		for(int j=0;j<num_rays;j++){
			int allinbox = 1;
			do{
				_a[3*j+0] = tools_uniform(-boxsize_x/2.,boxsize_x/2.);
				_a[3*j+1] = tools_uniform(-boxsize_x/2.,boxsize_x/2.);
				// Make sure line in completely in box.
				if (fabs(a[0]+b[0])>boxsize_x/2.) allinbox=0;
				if (fabs(a[0]-b[0])>boxsize_x/2.) allinbox=0;
				if (fabs(a[1]+b[1])>boxsize_y/2.) allinbox=0;
				if (fabs(a[1]-b[1])>boxsize_y/2.) allinbox=0;
			}while(allinbox=0);
			_a[3*j+2] = 0;
		}
#ifdef MPI
	}
	MPI_Bcast(_a,3*num_rays,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif //MPI

	for(int j=0;j<num_rays;j++){
		a[0] = _a[3*j+0];
		a[1] = _a[3*j+1];
		a[2] = _a[3*j+2];
		int transparency = 0;
		for (int i=0;i<root_n;i++){
#ifdef MPI
			int root_n_per_node = root_n/mpi_num;
			int proc_id = i/root_n_per_node;
			if (proc_id!=mpi_id) continue;
#endif //MPI
			transparency += tree_get_transparency(tree_root[i]);	
		}
		if (transparency>0) _t[j]++;
	}
	
#ifdef MPI
	int*   _tmax = malloc(num_rays*sizeof(int));
	MPI_Reduce(_t,_tmax,num_rays,MPI_INT,MPI_MAX,0,MPI_COMM_WORLD);
	free(_t);
	_t = _tmax;
	if (mpi_id==0){
#endif //MPI
		int num_intersecting = 0;
		for(int j=0;j<num_rays;j++){
			if (_t[j]>0) num_intersecting++;
		}
		FILE* of = fopen("transparency.txt","a+"); 
		if (of==NULL){
			printf("\n\nError while opening file.\n");
			return;
		}
		fprintf(of,"%e\t%e\t%e\t%0.8f\n",t,B,phi,(double)num_intersecting/(double)num_rays);
		fclose(of);
#ifdef MPI
	}
#endif //MPI
	free(_a);
	free(_t);
}

void problem_output(){
	for (double phi = 0; phi<M_PI; phi+=M_PI/10.){
		output_transparency(M_PI/2.,phi); //Normal to ring plane
		for (double B = M_PI/4.; B<M_PI/2.; B+=M_PI/20.){
			output_transparency(B,phi);
		}

	}
	if (output_check(10.*dt)){
		output_timing();
	}
	if (output_check(2.*M_PI/OMEGA)){
		char filename[256];
		sprintf(filename,"position_%08d.txt",position_id);
		output_ascii_mod(filename);
		sprintf(filename,"binary_%08d.txt",position_id);
		output_binary(filename);
		position_id++;
	}
}

void problem_finish(){
#ifdef MPI
	if (mpi_id==0){
#endif // MPI
		struct timeval tim;
		gettimeofday(&tim, NULL);
		double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
		output_double("runtime [s]",timing_final-timing_initial);
		system("cat config.log");
#ifdef MPI
	}
#endif // MPI
}
