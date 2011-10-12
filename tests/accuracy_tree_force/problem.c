#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "tree.h"
#include "particle.h"
#include "boundaries.h"
#include "output.h"
#include "gravity.h"


double energy_initial;

#ifdef GRAVITY_TREE
extern double opening_angle2;
#else // GRAVITY_TREE
double opening_angle2;
#endif // GRAVITY_TREE

void problem_init(int argc, char* argv[]){
	// Setup particle structures
	if (argc>1){
		opening_angle2 = atof(argv[1]);
		opening_angle2 *= opening_angle2;
	}
	boxsize 	= 2;
	dt 		= 1e-4;
	tmax 		= 10.*dt;
	int _N 		= 1000;
	init_box();
	// Initial conditions
	for (int i =0;i<_N;i++){
		struct particle p;
		double r;
		double rmax = 0.1*boxsize;
		do{
			p.x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
			p.y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
			p.z = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_z;
			r = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
		}while(r>rmax);
		double M = pow(r/rmax,3.);
		double vkep = sqrt(G*M/r);
		double phi = 2.*M_PI*((double)rand()/(double)RAND_MAX-0.5);
		double psi = 2.*M_PI*((double)rand()/(double)RAND_MAX-0.5);
		p.vx = vkep*sin(phi)*sin(psi)   +vkep*0.1*((double)rand()/(double)RAND_MAX-0.5);
		p.vy = vkep*cos(phi)*sin(psi)   +vkep*0.1*((double)rand()/(double)RAND_MAX-0.5);
		p.vz = vkep*cos(psi)            +vkep*0.1*((double)rand()/(double)RAND_MAX-0.5);
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m = 1./(double)_N;
		particles_add(p);
	}
	// No ghost boxes 
	nghostx = 0;
	nghosty = 0;
	nghostz = 0;
}
double average_relative_force_error(){	
	double error=0;
	double acceltot=0;
	for (int i=0; i<N; i++){
		double ax = particles[i].ax; 
		double ay = particles[i].ay; 
		double az = particles[i].az; 
		particles[i].ax = 0; 
		particles[i].ay = 0; 
		particles[i].az = 0; 
		for (int gbx=-nghostx; gbx<=nghostx; gbx++){
		for (int gby=-nghosty; gby<=nghosty; gby++){
		for (int gbz=-nghostz; gbz<=nghostz; gbz++){
			struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,gbz);
			for (int j=0; j<N; j++){
				if (i==j) continue;
				double dx = (gb.shiftx+particles[i].x) - particles[j].x;
				double dy = (gb.shifty+particles[i].y) - particles[j].y;
				double dz = (gb.shiftz+particles[i].z) - particles[j].z;
				double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
				double prefact = -G/(r*r*r)*particles[j].m;
				particles[i].ax += prefact*dx; 
				particles[i].ay += prefact*dy; 
				particles[i].az += prefact*dz; 
			}
		}
		}
		}
		acceltot += fabs((particles[i].ax));
		acceltot += fabs((particles[i].ay));
		acceltot += fabs((particles[i].az));
		error += fabs((particles[i].ax - ax));
		error += fabs((particles[i].ay - ay));
		error += fabs((particles[i].az - az));
	}
	printf("error = %f\n",acceltot);
	return error/acceltot;
}

void problem_inloop(){
}

void problem_output(){
}

void problem_finish(){
	FILE* of = fopen("error.txt","a+"); 
	tree_update();
	tree_update_gravity_data();
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_initloc = tim.tv_sec+(tim.tv_usec/1000000.0);
	for (int i=0;i<10;i++){
		gravity_calculate_acceleration();
	}
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error= average_relative_force_error();
	// These errors limits are empirical 
#ifdef QUADRUPOLE
	double error_limit = 0.03*pow(opening_angle2,4./2.);
#else  // QUADRUPOLE
	double error_limit = 0.054*pow(opening_angle2,3./2.);
#endif // QUADRUPOLE
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initloc);
#ifdef QUADRUPOLE
	fprintf(of,"Quadrupole, N = %d, opening_angle = %f",N,sqrt(opening_angle2));
#else
	fprintf(of,"Monopole, N = %d, opening_angle = %f",N,sqrt(opening_angle2));
#endif
	fprintf(of,"\n");
	fclose(of);
}
