#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
void collision_resolve_merge(struct reb_simulation* const mini, struct reb_collision c);
double E0;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();

    double planetesimal_mass = 1e-8;
    double amin = 0.45, amax = 0.75;        //for planetesimal disk
    double powerlaw = 0.5;
    
    int N_planetesimals = 110;
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HYBARID;
	//r->integrator	= REB_INTEGRATOR_IAS15;
	//r->integrator	= REB_INTEGRATOR_WHFAST;
    r->ri_hybarid.switch_ratio = 2;  //Hill radii
    r->ri_hybarid.CE_radius = 5.;          //X*radius
    r->testparticle_type = 1;
	r->heartbeat	= heartbeat;
    //r->usleep = 500;
    r->dt = 0.0001;

    r->collision = REB_COLLISION_DIRECT;
    //r->collision_resolve = collision_resolve_merge;
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(140);
    
    //planet 1
    {
        double a=0.5, m=5e-5, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p1.r = 1.6e-4;              //radius of particle is in AU!
        p1.id = r->N;
        reb_add(r, p1);
    }
    
    //planet 2
    {
        double a=0.7, m=5e-5, e=0.01, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = 1.6e-4;
        p2.id = r->N;
        reb_add(r, p2);
    }
    
    r->N_active = r->N;
    reb_move_to_com(r);
    
    //planetesimals
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(amin,amax,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.0001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, planetesimal_mass, a, 0., inc, Omega, apsis,phi);
		pt.r 		= 4e-5;
        pt.id = r->N;
		reb_add(r, pt);
    }
    
    system("rm -f energy.txt");
    E0 = reb_tools_energy(r);
    reb_integrate(r, INFINITY);
    
}

double tout = 1.;
void heartbeat(struct reb_simulation* r){
	if (reb_output_check(r, 100.*r->dt)){
		reb_output_timing(r, 0);
	}
    if (tout <r->t){
        tout *=1.01;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0+r->ri_hybarid.dE_offset)/E0);
        FILE* f = fopen("energy.txt","a+");
        int N_mini = 0;
        if (r->ri_hybarid.mini_active){
            N_mini = r->ri_hybarid.mini->N;
        }
        fprintf(f,"%e %e %d\n",r->t,relE,N_mini);
        fclose(f);
    }
}

//check for collisions in mini each heartbeat
void collision_resolve_merge(struct reb_simulation* const mini, struct reb_collision c){
    struct reb_simulation* r = mini->ri_hybarid.global;
    struct reb_particle* particles = mini->particles;
    int N_active = mini->N_active;
    int i = c.p1;
    int j = c.p2;
    
    if (i<N_active && j>=N_active){
        // Already ordered
    }else if (i>=N_active && j<N_active){
        // Will be done later
        return;
    }else{
        // Ignore collision between small particles
        return;
    }

    struct reb_particle* pi = &(particles[i]);
    struct reb_particle* pj = &(particles[j]);
                
    double invmass = 1.0/(pi->m + pj->m);
    double Ei = reb_tools_energy(mini);
    
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m += pj->m;
    
    reb_remove(mini,j,1);  //remove mini
    
    double Ef = reb_tools_energy(mini);
    r->ri_hybarid.dE_offset += Ei - Ef;
    printf("\n\tParticle %d collided with body %d from system at t=%f.\n",i,r->ri_hybarid.encounter_index[j],r->t);
    
    //remove from global and update global arrays
    int globalj = r->ri_hybarid.encounter_index[j];
    reb_remove(r,globalj,1);
    
    for(int k=globalj;k<r->N;k++){
        r->ri_hybarid.particles_prev[k] = r->ri_hybarid.particles_prev[k+1];
        r->ri_hybarid.is_in_mini[k] = r->ri_hybarid.is_in_mini[k+1];
    }
    r->ri_hybarid.encounter_index_N--;
    for(int k=j;k<r->ri_hybarid.encounter_index_N;k++) r->ri_hybarid.encounter_index[k] = r->ri_hybarid.encounter_index[k+1];
    for(int k=N_active;k<r->ri_hybarid.encounter_index_N;k++){
        if(r->ri_hybarid.encounter_index[k] > globalj) r->ri_hybarid.encounter_index[k]--; //1 fewer particles in index now
    }
}
