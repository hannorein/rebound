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
double lambda_crit;

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
	lambda_crit	 	= 4.*M_PI*M_PI*G*surface_density/OMEGA/OMEGA;
	boxsize 		= 12.*lambda_crit;
	int _N 			= (int)round(tau*boxsize*boxsize/(M_PI*particle_radius*particle_radius));
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.01*OMEGA*particle_radius;
	dt 			= 1e-2*2.*M_PI;
	tmax			= 20.*2.*M_PI;
	softening 		= 0.1*particle_radius;
	// Nothing needs to be changes below.
	init_box();
	printf("lambda_crit = %f\nsurface_density = %f\nboxsize = %f\nparticle_mass = %f\n",lambda_crit,surface_density,boxsize,particle_mass);
	// Initial conditions
	while(N<_N){
		struct particle p;
		p.x = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_x;
		p.y = ((double)rand()/(double)RAND_MAX-0.5)*boxsize_y;
		p.z = 3.*particle_radius*((double)rand()/(double)RAND_MAX-0.5);
		p.vx = 0;
		p.vy = -1.5*p.x*OMEGA;
		p.vz = 0;
		p.ax = 0;
		p.ay = 0;
		p.az = 0;
		p.m = particle_mass;
		p.r = particle_radius;
		particles_add(p);
	}
	// Do use ghost boxes in x and y
	nghostx = 3;
	nghosty = 3;
	nghostz = 0;
}

extern double collisions_plog;
double collisions_plog_last =0;
void calculate_viscosity(double* nutrans, double* nucoll, double* nugrav){
	double _nutrans =0;
	double _nugrav  =0;
	double mtot = 0;
	double sx2 = 3.*lambda_crit/2.;
	double sy2 = 3.*lambda_crit/2.;
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		mtot += p.m;
		double vyr = p.vy + 1.5 * OMEGA * p.x;
		_nutrans += p.m*p.vx*vyr;
		for (int gbx=0; gbx<=1; gbx++){
			for (int gby=-2; gby<=2; gby++){
				struct ghostbox gb = boundaries_get_ghostbox(gbx,gby,0);
				for (int j=0;j<N;j++){
					struct particle p2 = particles[j];
					if (p2.x<=p.x) continue;
					if (p2.x-p.x>sx2) continue;
					if (p.x-p2.x>sx2) continue;
					if (p2.y-p.y>sy2) continue;
					if (p.y-p2.y>sy2) continue;
					double dx = p2.x+gb.shiftx - p.x;
					double dy = p2.y+gb.shifty - p.y;
					double dz = p2.z - p.z;
					double r = sqrt(dx*dx + dy*dy + dz*dz);
					_nugrav += p.m*p2.m*dy*dx/(r*r*r);
				}
			}
		}
	}
	*nutrans = 2./(3.*OMEGA*mtot)*_nutrans;
	*nugrav  = 2./(3.*OMEGA*mtot)*_nugrav * (-G);
	*nucoll  = 2./(3.*OMEGA*mtot)*1./(t-collisions_plog_last)*collisions_plog;
}

double nutrans_tot=0, nucoll_tot=0, nugrav_tot=0; int nuN=0;

int reset = 1;
void problem_inloop(){
}
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
	output_timing();
	if (t>10.*2.*M_PI&&reset==1){
		collisions_plog_last 	= t;
		collisions_plog		= 0;
		nuN 			= 0;
		nutrans_tot 		= 0;
		nucoll_tot  		= 0; 
		nugrav_tot  		= 0;
		reset 			= 0;
	}
	if (output_check(1e-1*2.*M_PI/OMEGA)&&t>10.*2.*M_PI){
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
