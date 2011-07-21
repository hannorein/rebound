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
	boxsize = 2;
	softening = 0;//boxsize/100.;
	init_particles(10000);
	dt = 1e-4;
	tmax = dt/10.;
	// Initial conditions
	for (int i =0;i<N;i++){
		double r;
		double rmax = 0.1*boxsize;
		do{
			particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
			particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
			particles[i].z = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_z;
			r = sqrt(particles[i].x*particles[i].x+particles[i].y*particles[i].y+particles[i].z*particles[i].z);
		}while(r>rmax);
		double M = pow(r/rmax,3.);
		double vkep = sqrt(G*M/r);
		double phi = 2.*M_PI*((double)rand()/(double)RAND_MAX-0.5);
		double psi = 2.*M_PI*((double)rand()/(double)RAND_MAX-0.5);
		particles[i].vx = vkep*sin(phi)*sin(psi)   +vkep*0.1*((double)rand()/(double)RAND_MAX-0.5);
		particles[i].vy = vkep*cos(phi)*sin(psi)   +vkep*0.1*((double)rand()/(double)RAND_MAX-0.5);
		particles[i].vz = vkep*cos(psi)            +vkep*0.1*((double)rand()/(double)RAND_MAX-0.5);
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = 1./(double)N;
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
			struct ghostbox gb = get_ghostbox(gbx,gby,gbz);
			for (int j=N_active_first; j<N_active_last; j++){
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
	return error/acceltot;
}

void problem_inloop(){
}

void problem_output(){
}

void problem_finish(){
	FILE* of = fopen("error.txt","a+"); 
	calculate_forces();
	double error= average_relative_force_error();
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 1e-4;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
#ifdef QUADRUPOLE
	fprintf(of,"Quadrupole, N = %d, opening_angle = %f",N,sqrt(opening_angle2));
#else
	fprintf(of,"Monopole, N = %d, opening_angle = %f",N,sqrt(opening_angle2));
#endif
	fprintf(of,"\n");
	fclose(of);
}
