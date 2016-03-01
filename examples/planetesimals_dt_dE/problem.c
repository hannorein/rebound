/**
 * A.S. This is my planetesimal disk with special integrator for close encounters.
 *      Particle id's: 0 = star, 1 = massive body, 2 = planetesimal, 3 = CLOSE ENCOUNTER
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double E0, t_output, t_log_output;
char* mercury_dir; char* swifter_dir;

//temp
int output_xyz = 0, numdt=1, xyz_counter=0; double xyz_t = 36.76;//switch to 0 for no outputs
time_t t_ini;
int N_prev;
char* argv4;
char output_name[100] = {0};

//swifter/mercury compare
void output_to_mercury_swifter(struct reb_simulation* r, double HSR, double tmax, int n_output);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    int N_planetesimals = atoi(argv[2]);
    int seed = atoi(argv[3]);
    strcat(output_name,argv[4]); strcat(output_name,".txt"); argv4=argv[4];
    
	//Simulation Setup
	r->integrator	= REB_INTEGRATOR_HYBARID;
    r->ri_hybarid.CE_radius = 20.;         //X*radius
    r->testparticle_type = 1;
    r->heartbeat	= heartbeat;
    
    //which test?
    int test_dt = 0;
    double tmax;
    if(test_dt){//testing dt
        tmax = 1000 * 2 * 6.28319;   //1000 orbits of outer planet
        r->ri_hybarid.switch_ratio = 2;        //units of Hill radii
        r->dt = atof(argv[1]) * 6.28319;
    } else {//testing HSR
        tmax = 1e6 * 2 * 6.28319;   //1Myr orbits of outer planet
        r->ri_hybarid.switch_ratio = atof(argv[1]);
        r->dt = 0.01 * 6.28319; //100 dt's per orbital period
    }
    
    //collision
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->collisions_track_dE = 1;     //switch to track the energy from collisions/ejections
    
    //boundary
    //r->boundary			= REB_BOUNDARY_OPEN;
    //double boxsize = 5;
    //reb_configure_box(r, boxsize, 1, 1, 1);
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
    star.r		= 0.005;        // Radius of particle is in AU!
	reb_add(r, star);
   
    srand(seed);
    int n_output = 25000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = r->dt;
    printf("tlogoutput=%f",t_log_output);
    
    //planet 1
    {
        double a=1, m=5e-5, e=0, inc = reb_random_normal(0.00001);
        struct reb_particle p1 = {0};
        p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p1.r = 1.6e-4;              //radius of particle is in AU!
        p1.id = r->N;
        reb_add(r, p1);
    }
    
    //planet 2
    {
        double a=1.5, m=5e-5, e=0.01, inc=reb_random_normal(0.00001);
        struct reb_particle p2 = {0};
        p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p2.r = 1.6e-4;
        p2.id = r->N;
        reb_add(r, p2);
    }
    
    r->N_active = r->N;
    reb_move_to_com(r);
    
    //planetesimals
    double planetesimal_mass = 1e-8;
    double amin = 0.95, amax = 1.65;        //for planetesimal disk
    double powerlaw = 1;
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
    
    //energy
    E0 = reb_tools_energy(r);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,argv[4]); strcat(syss,"*");
    system(syss);
    
    t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    N_prev = r->N;
    argv4= argv[4];
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //elapsed time stuff
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    char timeout[200] = {0};
    strcat(timeout,argv[4]); strcat(timeout,"_elapsedtime"); strcat(timeout,".txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    
}

void heartbeat(struct reb_simulation* r){
    if(r->t > t_output){//log output
        t_output = r->t*t_log_output;
        
        double E = reb_tools_energy(r) + r->collisions_dE;
        double dE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("    dE=%e",dE);
        
        time_t t_curr = time(NULL);
        struct tm *tmp2 = gmtime(&t_curr);
        double time = t_curr - t_ini;
        
        FILE *append;
        append = fopen(output_name, "a");
        fprintf(append, "%.16f,%.16f,%d,%d,%.1f\n",r->t,dE,r->N,r->ri_hybarid.mini->N,time);
        fclose(append);
    }
    
    //record collisions in mini
    if(r->N < N_prev){
        char removed[200] = {0}; strcat(removed,argv4); strcat(removed,"_removed"); strcat(removed,".txt");
        FILE* append = fopen(removed,"a");
        fprintf(append,"Collision,%.5f\n",r->t);
        fclose(append);
        
        N_prev = r->N;
    }
    
    //ejections
    const double ED2 = 9;
    struct reb_particle* global = r->particles;
    struct reb_particle p0 = global[0];
    for(int i=1;i<r->N;i++){
        const double dx = global[i].x - p0.x;
        const double dy = global[i].y - p0.y;
        const double dz = global[i].z - p0.z;
        if(dx*dx+dy*dy+dz*dz > ED2){
            const double Ei = reb_tools_energy(r);
            reb_remove(r,i,1);
            const double Ef = reb_tools_energy(r);
            r->collisions_dE += Ei - Ef;
            
            char removed[200] = {0}; strcat(removed,argv4); strcat(removed,"_removed"); strcat(removed,".txt");
            FILE* append = fopen(removed,"a");
            fprintf(append,"Ejection,%.5f\n",r->t);
            fclose(append);
            
            N_prev = r->N;
        }
    }
    
    //xyz movie plots
    if(r->t >= xyz_t - 0.01*r->dt && output_xyz){
        double E = reb_tools_energy(r) + r->collisions_dE;
        double dE = fabs((E-E0)/E0);
        xyz_t += numdt*r->dt;
        //printf("\n%f,%d,%f,%f\n",xyz_t,numdt,r->dt,r->t);
        xyz_counter++;
        char filename[100]={0}; char ii[5] = {0};
        sprintf(ii, "%d", xyz_counter);
        strcat(filename,"xyz_outputs/hybarid"); strcat(filename,ii); strcat(filename,".txt");
        FILE* output = fopen(filename,"w");
        for(int i=0; i<r->N;i++) fprintf(output,"%f,%d,%.16f,%.16f,%.16f,%.16f\n",r->t,r->particles[i].id,r->particles[i].x,r->particles[i].y,r->particles[i].z,dE);
        fclose(output);
    }
}

