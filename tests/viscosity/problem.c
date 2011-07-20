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
double particle_radius = 1;
double r_hstar 		= 0.5;

void problem_init(int argc, char* argv[]){
	// Setup constants
	OMEGA 	= 1.;
	G	= 1;
	// Setup particle structures
	if (argc>1){
		r_hstar = atof(argv[1]);
	}
	double r_h 		= r_hstar*2.*particle_radius;
	double particle_mass 	= r_h*r_h*r_h*3./2.;
	double tau 		= 0.5;
	double surface_density 	= tau/(M_PI*particle_radius*particle_radius)*particle_mass;
	double lambda_crit 	= 4.*M_PI*M_PI*G*surface_density/OMEGA/OMEGA;
	boxsize 		= 150;
	N 			= (int)round(tau*boxsize*boxsize/(M_PI*particle_radius*particle_radius));
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.01*OMEGA*particle_radius;
	dt 			= 1e-2;
	tmax			= 20.*2.*M_PI;
	softening 		= 0.1*particle_radius;
	// Nothing needs to be changes below.
	init_particles(N);
	printf("lambda_crit = %f\nsurface_density = %f\nboxsize = %f\nparticle_mass = %f\n",lambda_crit,surface_density,boxsize,particle_mass);
	// Initial conditions
	for (int i =0;i<N;i++){
		particles[i].x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		particles[i].y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		particles[i].z = 3.*particle_radius*((double)rand()/(double)RAND_MAX-0.5);
		particles[i].vx = 0;
		particles[i].vy = -1.5*particles[i].x*OMEGA;
		particles[i].vz = 0;
		particles[i].ax = 0;
		particles[i].ay = 0;
		particles[i].az = 0;
		particles[i].m = particle_mass;
		particles[i].r = particle_radius;
	}
	// Do use ghost boxes in x and y
	nghostx = 3;
	nghosty = 3;
	nghostz = 0;
}

void problem_inloop(){

}

extern double collisions_plog;
double collisions_plog_last =0;
void calculate_viscosity(double* nutrans, double* nucoll, double* nugrav){
	double _nutrans =0;
	double _nugrav  =0;
	double mtot =0;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		mtot += p.m;
		double vyr = p.vy + 1.5 * OMEGA * p.x;
		_nutrans += p.m*p.vx*vyr;
		for (int j=0;j<N;j++){
			struct particle p2 = particles[j];
			if (p2.x<=p.x) continue;
			for (int gbx=0; gbx<=nghostx; gbx++){
				for (int gby=-nghosty; gby<=nghosty; gby++){
					struct ghostbox gb = get_ghostbox(gbx,gby,0);
					double dx = p2.x+gb.shiftx - p.x;
					double dy = p2.y+gb.shifty - p.y;
					double dz = p2.z - p.z;
					double r2 = dx*dx + dy*dy + dz*dz + softening*softening;
					_nugrav += p.m*p2.m*dy*dx/pow(r2,3./2.);
				}
			}
		}
	}
	*nutrans = 2./(3.*OMEGA*mtot)*_nutrans;
	*nugrav  = 2./(3.*OMEGA*mtot)*_nugrav * (-G);
	*nucoll  = 2./(3.*OMEGA*mtot)*1./(t-collisions_plog_last)*collisions_plog;
}

double nutrans_tot=0, nucoll_tot=0, nugrav_tot=0; int nuN=0;
void calculate_viscosity_tot(){
	double nutrans,nucoll,nugrav;
	calculate_viscosity(&nutrans,&nucoll,&nugrav);
	nutrans_tot +=nutrans;
	nucoll_tot  +=nucoll;
	nugrav_tot  +=nugrav;
	collisions_plog_last=t;
	collisions_plog = 0;
	nuN++;
}


void problem_output(){
	if (output_check(1e-1*2.*M_PI/OMEGA)&&t>0){
		calculate_viscosity_tot();
	}
}


void problem_finish(){
	calculate_viscosity_tot();
	FILE* of = fopen("error.txt","a+"); 
	double error= 1;  // TODO!
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing_final = tim.tv_sec+(tim.tv_usec/1000000.0);
	double error_limit = 0.1;
	int pass = (error>error_limit?0:1);
	fprintf(of,"%d\t%e\t%e\t%e\t%d\t%e\t",N,dt,error,error_limit,pass,timing_final-timing_initial);
	fprintf(of,"N_init = %d",0);
	fprintf(of,"\n");
	fclose(of);
	
	FILE* ofv = fopen("viscosity.txt","a"); 
	fprintf(ofv,"%e\t%e\t%e\t%e\t%e\n",t,r_hstar,nutrans_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN,nucoll_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN,nugrav_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN);
	fclose(ofv);

}
