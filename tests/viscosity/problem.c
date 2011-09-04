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

extern double OMEGA;
extern double coefficient_of_restitution;
extern double minimum_collision_velocity;
double particle_radius = 1;
double r_hstar 		= 0.5;
double lambda_crit;
double sx2;
double sy2;

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
	sx2 			= 3.*lambda_crit/2.;
	sy2 			= 3.*lambda_crit/2.;
	if (argc>2){
		boxsize 	=lambda_crit*atof(argv[2]);
	}
	int _N 			= (int)round(tau*boxsize*boxsize/(M_PI*particle_radius*particle_radius));
	coefficient_of_restitution = 0.5;
	minimum_collision_velocity = 0.01*OMEGA*particle_radius;
	dt 			= 1e-2*2.*M_PI;
	tmax			= 40.*2.*M_PI;
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


double calculate_viscosity_trans(){
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
		
double calculate_viscosity_grav_for_particle_from_cell(int i, struct ghostbox gb,struct cell* node){
	double nu = 0;
	double dx = gb.shiftx - node->mx;
	double dy = gb.shifty - node->my;
	double dz = gb.shiftz - node->mz;
	if (node->pt<0){
		if (fabs(dx)-node->w/2. < sx2 && fabs(dy)-node->w/2. < sy2){
			for(int o=0;o<8;o++){
				if (node->oct[o]!=NULL){
					nu += calculate_viscosity_grav_for_particle_from_cell(i,gb,node->oct[o]);
				}
			}
		} 
	}else{
		if (fabs(dx) < sx2 && fabs(dy) < sy2 && dx > 0 && i!=node->pt){
			double r = sqrt(dx*dx + dy*dy + dz*dz);
			nu += particles[i].m*particles[node->pt].m*dy*dx/(r*r*r);
		}
	}
	return nu;
}

double calculate_viscosity_grav_for_particle(int pt, struct ghostbox gb) {
	double nu =0;
	for(int i=0;i<root_n;i++){
		struct cell* node = tree_root[i];
		if (node!=NULL){
			nu += calculate_viscosity_grav_for_particle_from_cell(pt, gb, node);
		}
	}
	return nu;
}

double calculate_viscosity_grav(){
	double nu  	= 0;
	double mtot 	= 0;
	for (int i=0;i<N;i++){
		mtot	+= particles[i].m;
	}

	// Summing over all Ghost Boxes
	for (int gbx=-nghostx; gbx<=nghostx; gbx++){
	for (int gby=-nghosty; gby<=nghosty; gby++){
		struct ghostbox gb_cache = boundaries_get_ghostbox(gbx,gby,0);
		// Summing over all particle pairs
		for (int i=0; i<N; i++){
			struct ghostbox gb = gb_cache;
			// Precalculated shifted position
			gb.shiftx += particles[i].x;
			gb.shifty += particles[i].y;
			gb.shiftz += particles[i].z;
			nu += calculate_viscosity_grav_for_particle(i, gb);
		}
	}}
	return 2./(3.*OMEGA*mtot)*nu * (-G);
}

double nutrans_tot=0, nucoll_tot=0, nugrav_tot=0; int nuN=0;
int reset = 1;

void calculate_viscosity(){
	nutrans_tot 	+= calculate_viscosity_trans();
	nugrav_tot 	+= calculate_viscosity_grav();
	nucoll_tot	 = calculate_viscosity_coll();
	nuN++;
}

void problem_inloop(){
}

void problem_output(){
	output_timing();
	if (t>tmax/2.){
      		if(reset==1){
			collisions_plog_last 	= t;
			collisions_plog		= 0;
			nuN 			= 0;
			nutrans_tot 		= 0;
			nucoll_tot  		= 0; 
			nugrav_tot  		= 0;
			reset 			= 0;
		}else{
			calculate_viscosity();
		}
	}
}


void problem_finish(){
	FILE* ofv = fopen("viscosity.txt","a"); 
	double nu_trans = nutrans_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN;
	double nu_coll 	= nucoll_tot/(particle_radius*particle_radius*OMEGA);
	double nu_grav	= nugrav_tot/(particle_radius*particle_radius*OMEGA)/(double)nuN;
	fprintf(ofv,"%e\t%e\t%e\t%e\t%e\n",t,r_hstar,nu_trans,nu_coll,nu_grav);
	fclose(ofv);

}
