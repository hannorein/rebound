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
#include "integrator_whfast.h"
#include "tools.h"
#include "../examples/planetesimals2/functions.h"

void heartbeat(struct reb_simulation* r);
char plntdir[200] = "output/planet_", lgnddir[200] = "output/planet_", removeddir[200]="output/planet_", CEprint[200]="output/planet_";

double tmax, planetesimal_mass, E0, n_output, dt_ini, t_output, t_log_output, ias_timestep, soft, dE_collision = 0, movie_ti, movie_tf, ejection_distance2, dE_prev=0;
int N_encounters = 0, N_encounters_previous, N_encounters_tot = 0, HYBRID_ON, n_o=0, movie_mini = 0, movie_output, movie_counter, movie_mc = 0, movie_output_interval, movie_output_file_per_time;
int* encounter_index; int* previous_encounter_index; double* Hill2; double* x_prev; double* y_prev; double* z_prev; double t_prev;
struct reb_simulation* s; struct reb_simulation* r;

int main(int argc, char* argv[]){
    //System constants
    tmax = atof(argv[1]);
    int N_planetesimals = atoi(argv[2]);
    //double M_planetesimals = 3e-6;                        //Tot. Mass of all planetesimals (Earth mass, 3e-6)
    //planetesimal_mass = M_planetesimals/N_planetesimals;  //mass of each planetesimal
    planetesimal_mass = 3e-8;                               //each is a moon
    double M_planetesimals = planetesimal_mass*N_planetesimals;
    double ias_epsilon = 1e-8;                              //sets precision of ias15
    double HSR2 = 3;                                        //Transition boundary bet. WHFAST & IAS15. Units of Hill^2
    double dRHill = 0.125;                                   //Sets the timestep - max # Hill radii/timestep.
    soft = 1.6e-4/10.;                                          //gravity softening length scale in AU. R_Neptune/100.
    int seed = atoi(argv[3]);
    
    //switches
    HYBRID_ON = 1;
    int turn_planetesimal_forces_on = 1;
    int p1_satellite_on = 0;
    int mercury_swifter_output = 1;
    //movie
    movie_output = 0;
    movie_output_interval = 1;                 //number of dt per movie output (used for swifter/mercury too)
    if(movie_output == 1){
        movie_mini = 1;                         //whether to output particles from mini or global
        movie_output_file_per_time = 1;          //output a file every unit of time (best for analyzing mini), or output file for each body and keep filling them
        movie_ti = 0;
        movie_tf = 1;
        system("rm -v movie/movie_output/*.txt");
        movie_counter = movie_output_interval;
    }
    
	//Simulation Setup
    r = reb_create_simulation();
	r->integrator	= 1;                                  //REB_INTEGRATOR_IAS15 = 0, WHFAST = 1, WH=3, HYBRID = 5
	r->collision	= REB_COLLISION_NONE;
	r->heartbeat	= heartbeat;
    r->ri_hybrid.switch_ratio = HSR2;
    r->softening = soft;
    if(turn_planetesimal_forces_on==1)r->additional_forces = planetesimal_forces_global;
    //r->ri_whfast.corrector 	= 11;
    //r->usleep   = 5000; //larger the number, slower OpenGL simulation
    
    //Boundary stuff
    //r->boundary     = REB_BOUNDARY_OPEN;
    //double boxsize = 20;
    //reb_configure_box(r, boxsize, 3, 3, 1);
    
    srand(seed);        //This needs to be here, just before rand() is called
    
	// Initial conditions
	struct reb_particle star = {0};
	star.m 		= 1;
	star.r		= 0.005;        // Radius of particle is in AU!
    star.id     = 0;            // 0 = star
	reb_add(r, star);
    
    double amin, amax;  //for planetesimal disk
    
    //planet 1
    double a=0.5, m=5e-4, e=0, inc = reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
    p1.r = 1.6e-4;              //radius of particle is in AU!
    p1.id = r->N;
    reb_add(r, p1);
    dt_ini = calc_dt(r, m, star.m, a, dRHill, 1);
    
    amin = a;
    
    //planet 2
    a=0.7, m=5e-4, e=0.01, inc=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
    p2.r = 1.6e-4;
    p2.id = r->N;
    reb_add(r, p2);
    dt_ini = calc_dt(r, m, star.m, a, dRHill, dt_ini);

    /*
    double a=5.2, m=0.0009543, e=0, inc = reb_random_normal(0.00001);
    struct reb_particle p1 = {0};
    p1 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0.85);
    p1.r = 0.00046732617;
    p1.id = r->N;
    reb_add(r, p1);
    dt_ini = calc_dt(r, m, star.m, a, dRHill, 1);
    
    amin = a;
    
    //planet 2
    a=9.5, m=0.0002857, e=0.0, inc=reb_random_normal(0.00001);
    struct reb_particle p2 = {0};
    p2 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
    p2.r = 0.000389256877;
    p2.id = r->N;
    reb_add(r, p2);
    dt_ini = calc_dt(r, m, star.m, a, dRHill, dt_ini);
     
    a=19.2, m=0.00004365, e=0.0, inc=reb_random_normal(0.00001);
    struct reb_particle p3 = {0};
    p3 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
    p3.r = 0.000169534499;
    p3.id = r->N;
    reb_add(r, p3);
    dt_ini = calc_dt(r, m, star.m, a, dRHill, dt_ini);
     
    a=30.1, m=5e-5, e=0.0, inc=reb_random_normal(0.00001);
    struct reb_particle p4 = {0};
    p4 = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
    p4.r = 1.6e-4;
    p4.id = r->N;
    reb_add(r, p4);
    dt_ini = calc_dt(r, m, star.m, a, dRHill, dt_ini);
    */
    amax = a;
    ejection_distance2 = pow(3*amax,2);     //distance at which particles removed from simulation (squared)
    
    //calc dt
    if(r->integrator == REB_INTEGRATOR_IAS15){
        dt_ini /= 3.;
        printf("dt = %f \n",dt_ini);
        dRHill = -1;
    }
    char timebuff[32] = {0};
    sprintf(timebuff, "%e", dt_ini);
    r->dt = atof(timebuff);
    printf("timesetep is dt = %.16f, ri_hybrid.switch_ratio=%f \n",r->dt,r->ri_hybrid.switch_ratio);
    
    //N_active and move to COM
    r->N_active = r->N;
    if(r->integrator != REB_INTEGRATOR_WH) reb_move_to_com(r);
    
    //Outputting points
    n_output = 10000;
    t_log_output = pow(tmax + 1, 1./(n_output - 1));
    t_output = r->dt;  //general output
    if(movie_output == 0)movie_output_interval = tmax/(r->dt*n_output);  //Number of dt per movie output (used for swifter/mercury too!!)
    
    //orbiting planetesimal/satellite
    if(p1_satellite_on == 1){
        double x=0.01;
        struct reb_particle pt = {0};
        //pt = reb_tools_orbit_to_particle(r->G, p1, 0, x, 0, 0, 0, 0, 0.1);    //works well with m2=5e-4
        //pt.y += p1.y;
        pt = reb_tools_orbit_to_particle(r->G, star, 0, a - x, 0, 0, 0, 0, -0.1); //m=planetesimal_mass?
        pt.r = 4e-5;            //I think radius of particle is in AU!
        pt.id = r->N;              //1 = planet
        reb_add(r, pt);
    }
    
    //planetesimals
    if(turn_planetesimal_forces_on == 0) planetesimal_mass = 0;
    double planetesimal_buffer = 0.1;   //max +/- distance from massive bodies, Chatterjee & Ford use 0.01
    double inner = amin - planetesimal_buffer, outer = amax + planetesimal_buffer, powerlaw = 0.5;
    while(r->N<N_planetesimals + r->N_active){
		struct reb_particle pt = {0};
		double a	= reb_random_powerlaw(inner,outer,powerlaw);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        double inc = reb_random_normal(0.0001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, 0, a, 0, inc, Omega, apsis,phi);
		pt.r 		= 4e-5;
        pt.id       = r->N;
		reb_add(r, pt);
	}
    
    //Ini malloc arrays
    x_prev = calloc(sizeof(double),r->N);           //Previous global positions for interpolating
    y_prev = calloc(sizeof(double),r->N);
    z_prev = calloc(sizeof(double),r->N);
    encounter_index = malloc(sizeof(int));          //encounter index
    previous_encounter_index = malloc(sizeof(int));
    Hill2 = calloc(sizeof(double),r->N);             //Hill radius squared for fast calc.
    calc_Hill2(r);
    
    //input files for swifter/mercury
    if(mercury_swifter_output == 1) output_to_mercury_swifter(r, sqrt(HSR2), tmax, n_output, movie_output_interval);
    
    //Initializing stuff
    legend(plntdir, lgnddir, removeddir, CEprint, r, tmax, planetesimal_mass, M_planetesimals, N_planetesimals,inner, outer, powerlaw, star.m, dRHill,ias_epsilon,seed,HYBRID_ON);
    E0 = calc_Etot(r, soft, 0);
    
    //Ini mini
    s = reb_create_simulation();    //initialize mini simulation (IAS15)
    double ias_subtime = 3.;        //sets how much smaller ias15 vs. global timestep is
    ias_timestep = r->dt/ias_subtime;
    ini_mini(r,s,ias_epsilon,turn_planetesimal_forces_on,ias_timestep,soft);
    clock_t t_ini = clock_start();
    
    //Integrate!
    reb_integrate(r, tmax);
    
    //finish
    clock_finish(t_ini,N_encounters_tot,r->N - r->N_active,lgnddir);
    printf("N_outputs total=%d\n",n_o);
    global_free();
}

void heartbeat(struct reb_simulation* r){
    double min_r = 1e8, max_val = 1e-8; int output_it = 0;
    if(HYBRID_ON == 1){
        if(N_encounters_previous == 0){
            check_for_encounter(r, s, &N_encounters, N_encounters_previous, &min_r, &max_val, removeddir, &output_it, &dE_collision, soft, ejection_distance2);
            if(N_encounters > 0){//1st update in a while, update mini massive bodies, add particles, no int
                s->t = r->t;
                int N_active = s->N_active;
                struct reb_particle* global = r->particles;
                struct reb_particle* mini = s->particles;
                for(int i=0; i<N_active; i++) mini[i] = global[i];
                add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,CEprint,soft,dE_collision,E0); //remove soft,dE_collision,E0 later
                update_previous_global_positions(r, N_encounters);
            } //otherwise do nothing.
        } else { //integrate existing mini, update global, add/remove new/old particles.
            reb_integrate(s, r->t);
            update_global(s,r,N_encounters_previous);
            check_for_encounter(r, s, &N_encounters, N_encounters_previous, &min_r, &max_val, removeddir, &output_it, &dE_collision, soft, ejection_distance2);
            add_or_subtract_particles(r,s,N_encounters,N_encounters_previous,CEprint,soft,dE_collision,E0); //remove soft,dE_collision,E0 later
            update_previous_global_positions(r, N_encounters);
        }
        update_encounter_indices(&N_encounters, &N_encounters_previous);
    }
    
    double E1 = calc_Etot(r, soft, dE_collision);
    double dE = fabs((E1 - E0)/E0);
    double comp = 0;
    if(r->t > r->dt) comp = dE/dE_prev;
    if(comp > 2 && r->t > 100){
        fprintf(stderr,"\n\033[1mEnergy Error Jumped by %fx \033[0m at t=%f\n",comp,r->t);

        FILE *append;
        append = fopen(plntdir, "a");
        fprintf(append, "%.16f,%.16f, %d, %.12f,%.12f,%.16f,-20000\n",r->t,s->t,r->N,min_r,max_val,dE);
        fclose(append);
    } else if(output_it == 1){
        FILE *append;
        append = fopen(plntdir, "a");
        fprintf(append, "%.16f,%.16f, %d, %.12f,%.12f,%.16f,-10000\n",r->t,s->t,r->N,min_r,max_val,dE);
        fclose(append);
    }
    dE_prev = dE;
    
    //OUTPUT stuff*******
    if(r->t > t_output || r->t <= r->dt){
        t_output = r->t*t_log_output;
        n_o++;
        FILE *append;
        append = fopen(plntdir, "a");
        fprintf(append, "%.16f,%.16f, %d, %.12f,%.12f,%.16f,%d\n",r->t,s->t,r->N,min_r,max_val,dE,n_o);
        fclose(append);
        
        reb_output_timing(r, 0);    //output only when outputting values. Saves some time
    }
    
    /*
    if(r->t < 0.05){
        struct reb_particle* g = r->particles;
        double dx = g[1].x - g[159].x;
        double dy = g[1].y - g[159].y;
        double dz = g[1].z - g[159].z;
        double rr = sqrt(dx*dx + dy*dy + dz*dz);
        printf("\n t=%f, r=%f\n",r->t,rr);
    }*/
    
    //output movie - outputs in heliocentric coords
    if(movie_output == 1){
        if(r->t > movie_ti && r->t < movie_tf){
            if(movie_counter >= movie_output_interval){
                char* dir = "movie/movie_output/hybridbody";
                struct reb_particle* particles; int N; double t;
                if(movie_mini == 1){particles=s->particles; N=s->N; t=s->t;} else {particles=r->particles; N=r->N; t=r->t;}
                if(movie_output_file_per_time == 1)output_frame_per_time(particles, dir, N, t, &movie_mc); else output_frame_per_body(particles, dir, N, t);
                movie_counter = 0;
            }
            movie_counter++;
        }
    }
    
}

/*if(N_encounters == 0){//last particle just leaving
 printf("\n");
 struct reb_particle comr = reb_get_com(r);
 struct reb_particle coms = reb_get_com(s);
 printf("com: dx=%.16f, dy=%.16f, dz=%.16f, dvx=%.16f, dvy=%.16f, dvz=%.16f\n",comr.x - coms.x, comr.y - coms.y,comr.z - coms.z,comr.vx - coms.vx,comr.vy - coms.vy,comr.vz - coms.vz);
 struct reb_particle* global = r->particles;
 struct reb_particle* mini = s->particles;
 for(int i=0;i<r->N_active;i++){
 printf("%d, dx=%.16f, dy=%.16f, dz=%.16f, dvx=%.16f, dvy=%.16f, dvz=%.16f\n",global[i].id, global[i].x-mini[i].x,global[i].y-mini[i].y,global[i].z-mini[i].z,global[i].vx-mini[i].vx,global[i].vy-mini[i].vy,global[i].vz-mini[i].vz);
 printf("%d, px=%.16f, py=%.16f, pz=%.16f, pvx=%.16f, pvy=%.16f, pvz=%.16f\n",global[i].id, (global[i].x-mini[i].x)/global[i].x,(global[i].y-mini[i].y)/global[i].y,(global[i].z-mini[i].z)/global[i].z,(global[i].vx-mini[i].vx)/global[i].vx,(global[i].vy-mini[i].vy)/global[i].vy,(global[i].vz-mini[i].vz)/global[i].vz);
 printf("%d, x=%.16f, y=%.16f, z=%.16f, vx=%.16f, vy=%.16f, vz=%.16f\n",global[i].id, global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz);
 }
 exit(0);
 }*/

/*
 //output error stuff - every iteration
 if(fabs((E1 - E0)/E0) > 1e-6){
 if(err_print_msg == 0){
 err_print_msg++;
 fprintf(stderr,"\n\033[1mERROR EXCEEDED for %s\033[0m, t=%.16f.\n",plntdir,r->t);
 }
 }*/

/*
 for(int i=0;i<r->N;i++){
 if(fabs(global[i].ax) + fabs(global[i].ay) + fabs(global[i].az) > 100){
 for(int j=0;j<N_encounters_previous;j++){
 if(global[i].id == previous_encounter_index[j]){
 printf("mini integrated large acc, par %d: ax=%f,ay=%f,az=%f\n",global[i].id,global[i].ax,global[i].ay,global[i].az);
 }
 }
 }
 }*/

/*
 if(reb_output_check(r,tmax/100)){
 FILE *xyz_output;
 xyz_output = fopen(removeddir, "a");
 struct reb_particle* global = r->particles;
 fprintf(xyz_output, "%.16f\n",r->t);
 for(int i=0;i<r->N;i++){
 fprintf(xyz_output, "%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
 }
 fclose(xyz_output);
 }*/

/*
 if(r->t > 74.5){
 FILE *xyz_output;
 xyz_output = fopen(removeddir, "a");
 struct reb_particle* global = r->particles;
 //fprintf(xyz_output, "%.16f,rmin=%.16f,vmax/rmin=%.16f\n",r->t,min_r,max_val);
 int i = 20;
 int j=1;
 //fprintf(xyz_output, "particle %d,x=%.16f,y=%.16f,z=%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
 //fprintf(xyz_output, "%.16f,%.16f,%.16f,%.16f,%.16f\n",r->t,s->t,fabs(global[i].x-dxold1),fabs(global[i].y-dyold1),fabs(global[i].z-dzold1));
 //fprintf(xyz_output, "planet %d,x=%.16f,y=%.16f,z=%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",j,global[j].x,global[j].y,global[j].z,global[j].vx,global[j].vy,global[j].vz,global[j].ax,global[j].ay,global[j].az);
 fprintf(xyz_output, "%.16f,%.16f,%.16f,%.16f,%.16f\n",r->t,s->t,fabs(global[j].x-dxold2),fabs(global[j].y-dyold2),fabs(global[j].z-dzold2));
 dxold1 = global[i].x; dyold1 = global[i].y; dzold1 = global[i].z;
 dxold2 = global[j].x; dyold2 = global[j].y; dzold2 = global[j].z;
 //for(int i=0;i<r->N;i++){
 //    fprintf(xyz_output, "%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f\n",i,global[i].x,global[i].y,global[i].z,global[i].vx,global[i].vy,global[i].vz,global[i].ax,global[i].ay,global[i].az);
 //}
 fclose(xyz_output);
 }*/


//OUTPUT stuff*******
/*double E_curr = 0, K_curr = 0, U_curr = 0, L_curr = 0, a_p = 0, d_p = 0, e_p = 0, t = r->t;
 calc_ELtot(&E_curr, &K_curr, &U_curr, &L_curr, planetesimal_mass, r); //calcs Etot all in one go.
 for(int i=1;i<r->N_active;i++){
 calc_ae(&a_p, &e_p, &d_p, r, i, t);
 
 FILE *append;
 append=fopen(plntdir, "a");
 fprintf(append,"%f,%.8f,%.8f,%.16f,%.16f,%.16f,%.16f,%.16f\n",t,a_p,e_p,fabs((E_ini - E_curr)/E_ini),fabs((K_ini - K_curr)/K_ini), fabs((U_ini - U_curr)/U_ini),fabs((L_ini - L_curr)/L_ini),d_p);
 fclose(append);
 E_curr = E_ini; L_curr = L_ini;
 }*/

#include "functions.c"
