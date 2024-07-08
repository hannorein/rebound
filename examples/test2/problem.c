/**
 * Dust evolution after DART impact
 *
 * This example shows how to integrate dust particles
 * using the IAS15 integrator with additional forces.
 * The example sets the function pointer `additional_forces`
 * to a function that describes the radiation forces and
 * non-gravitational perturbations of Dimorphos and Didymos.
 *
 * The output is custom too, as given in the function  `heartbeat`.
 * 
 * Written by Yun Zhang (2023 May)
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include <ejecta_transform.h>

void force_radiation(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

const double AU = 1.495978707e11;
//double Q_pr = 1.0;  // reflectivity coefficient of radiation pressure
const double tmax = 300*24*3600;  // 300 days
const double mass_system = 5.5e11; // kg
const double sep_system = 1170.0; // seperation m
const double vol_didy = 0.2295409644951028;  // km^3
const double vol_dimor = 0.001830603200702610;
const double r_Didy_Bary_x = 155002790587.400; /* 2022-Sep-27 23:14 (https://ssd.jpl.nasa.gov/horizons/app.html#/) */ //双小行星系统质心相对于太阳系黄道面的位置
const double r_Didy_Bary_y = 16396517401.5340;
const double r_Didy_Bary_z = -8556549730.13832;
const double v_Didy_Bary_x = -7790.84422295264;    //双小行星系统质心相对于太阳系黄道面的速度
const double v_Didy_Bary_y = 33143.3809500845;
const double v_Didy_Bary_z = 1022.84858862050;
const double T11 = -0.182453930731996;
const double T12 = 0.971278608867291;
const double T13 = -0.152736462959145;
const double T21 = 0.971278608867291;
const double T22 = 0.202182756110269;
const double T23 = 0.125459145097038;
const double T31 = 0.152736462959145;
const double T32 = -0.125459145097038;
const double T33 = -0.980271174621727;
const double r_dust = 0.001; // dust particle radius, m -> require to change SRP_coe as well!!!
const double rho_dust = 3000; // dust particle density, kg/m^3
const double Rsq_didy = 850.0/2.0 * 850.0/2.0;
const double Rsq_dimor = 175.0/2.0 * 175.0/2.0;
const double Rsq_hill = 70000.0*70000.0;  // twice Hill radius of D-D system, m

// J2
const double J2_didy = 0.0956324486828653;  // J2 of Didymos
const double J2_dimor = 0.113929814552540;  // J2 of Dimorphos

// SRP
const double c = 2.99792458e8;         // speed of light.
//double Fsun = 1367.0;  /* integrated stellar flux at 1 au, W/m^2 */
const double SRP_coe = 1.0 * 1367.0/2.99792458e8 * 3.0/4.0/3000; // Q_pr * Fsun/c * 3.0/4.0/rho_dust

int main(int argc, char* argv[]){
    
    // Set the number of OpenMP threads to be the number of processors
    // omp_set_num_threads(4);
    
    // Setup simulation structure and 3D visualization server
    struct reb_simulation* r = reb_simulation_create();
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->dt             = 10.0;    // Initial timestep, s
    r->N_active        = 3;     // Only the Sun and the Didymos-Dimorphos system are massive, the dust particles are treated as test particles
    r->additional_forces     = force_radiation;
    r->heartbeat         = heartbeat;
    r->G = 6.6743e-11; //  m^3 / kg s^2
    

    // Didymos
    double mass_didy = mass_system*vol_didy/(vol_didy+vol_dimor);
    double mass_dimor = mass_system*vol_dimor/(vol_didy+vol_dimor);
    struct reb_particle Didymos = {0};
    double r_didy_com = -vol_dimor*sep_system/(vol_didy+vol_dimor); // distance of Didymos to center of mass of the system
    double v_didy_com = sqrt(-r->G*mass_dimor/pow(sep_system,2.0)*r_didy_com);
    Didymos.m  = mass_didy;
    Didymos.r = 850.0/2.0;
    Didymos.x  = r_didy_com;
    Didymos.y  = 0.0;
    Didymos.z  = 0.0;
    Didymos.vx =              T11 * v_Didy_Bary_x + T12 * v_Didy_Bary_y + T13 * v_Didy_Bary_z;
    Didymos.vy = v_didy_com + T21 * v_Didy_Bary_x + T22 * v_Didy_Bary_y + T23 * v_Didy_Bary_z;
    Didymos.vz =              T31 * v_Didy_Bary_x + T32 * v_Didy_Bary_y + T33 * v_Didy_Bary_z;
    Didymos.hash = 1;
    reb_simulation_add(r, Didymos);
    
    // Dimorphos
    struct reb_particle Dimorphos = {0};
    double r_dimor_com = vol_didy*sep_system/(vol_didy+vol_dimor);
    double v_dimor_com = sqrt(r->G*mass_didy/pow(sep_system,2.0)*r_dimor_com);
    Dimorphos.m  = mass_dimor;
    Dimorphos.r = 175.0/2.0;
    Dimorphos.x  = r_dimor_com;
    Dimorphos.y  = 0.0;
    Dimorphos.z  = 0.0;
    Dimorphos.vx =               T11 * v_Didy_Bary_x + T12 * v_Didy_Bary_y + T13 * v_Didy_Bary_z;
    Dimorphos.vy = v_dimor_com + T21 * v_Didy_Bary_x + T22 * v_Didy_Bary_y + T23 * v_Didy_Bary_z;
    Dimorphos.vz =               T31 * v_Didy_Bary_x + T32 * v_Didy_Bary_y + T33 * v_Didy_Bary_z;
    Dimorphos.hash = 2;
    reb_simulation_add(r, Dimorphos);
    
    // Sun
    struct reb_particle star = {0};
    star.m  = 1.9884e30;
    star.x = - T11 * r_Didy_Bary_x - T12 * r_Didy_Bary_y - T13 * r_Didy_Bary_z;
    star.y = - T21 * r_Didy_Bary_x - T22 * r_Didy_Bary_y - T23 * r_Didy_Bary_z;
    star.z = - T31 * r_Didy_Bary_x - T32 * r_Didy_Bary_y - T33 * r_Didy_Bary_z;
    star.hash = 3;
    reb_simulation_add(r, star);

    // Dust particles
    unsigned int N_particles = 3; // current number of particles (didy, dimor, sun)
    FILE *dust_file = fopen("/nuke/linfel/Ejecta/data_0Pa_160s.txt", "r");
    if (dust_file == NULL) {
        fprintf(stderr, "Error: Could not open file.\n");
        return 1;
    }

    ReadParticle rp; ////////////////

    while (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&rp.ID, &rp.x, &rp.y, &rp.z, &rp.vx, &rp.vy, &rp.vz, &rp.mass, &rp.density) == 9) {
	transform(&rp, r_dimor_com); ///////////////////
	
	struct reb_particle p = {0};
	p.m = rp.mass;
	p.r = r_dust;  ///////////////modifying SRP_coe is also required?
	
	p.x = rp.x;
	p.y = rp.y;
	p.z = rp.z;
	p.vx = rp.vx + T11 * v_Didy_Bary_x + T12 * v_Didy_Bary_y + T13 * v_Didy_Bary_z;
	p.vy = rp.vy + T21 * v_Didy_Bary_x + T22 * v_Didy_Bary_y + T23 * v_Didy_Bary_z;
	p.vz = rp.vz + T31 * v_Didy_Bary_x + T32 * v_Didy_Bary_y + T33 * v_Didy_Bary_z;
	
	// check if the particle collides with Didymos or Dimorphos
	if ( ( (p.x-Didymos.x)*(p.x-Didymos.x) + (p.y-Didymos.y)*(p.y-Didymos.y) + (p.z-Didymos.z)*(p.z-Didymos.z) ) < Rsq_didy )
	    continue;
	if ( ( (p.x-Dimorphos.x)*(p.x-Dimorphos.x) + (p.y-Dimorphos.y)*(p.y-Dimorphos.y) + (p.z-Dimorphos.z)*(p.z-Dimorphos.z) ) < Rsq_dimor )
	    continue;
	
	N_particles++;
	p.hash = N_particles;
	
	reb_simulation_add(r, p);
    }
    fclose(dust_file);
 
    fprintf(stdout, "Total particle number: %i\n", N_particles);
    
    //reb_move_to_com(r);

    system("rm -v particles.txt");
    system("rm -v collide.txt");

    reb_simulation_integrate(r, tmax);
    fprintf(stdout, "\n");
}

void force_radiation(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    const struct reb_particle Didymos = particles[0];
    const struct reb_particle Dimorphos = particles[1];
    const struct reb_particle star = particles[2];            // cache
    const int N = r->N;
    
    // Sun-Didymos vector
    double prx_didy_star = Didymos.x-star.x;
    double pry_didy_star = Didymos.y-star.y;
    double prz_didy_star = Didymos.z-star.z;
    double pr = 1.0/sqrt(prx_didy_star*prx_didy_star + pry_didy_star*pry_didy_star + prz_didy_star*prz_didy_star);
    prx_didy_star *= pr;
    pry_didy_star *= pr;
    prz_didy_star *= pr;
    
    // Sun-Dimorphos vector
    double prx_dimor_star = Dimorphos.x-star.x;
    double pry_dimor_star = Dimorphos.y-star.y;
    double prz_dimor_star = Dimorphos.z-star.z;
    pr = 1.0/sqrt(prx_dimor_star*prx_dimor_star + pry_dimor_star*pry_dimor_star + prz_dimor_star*prz_dimor_star);
    prx_dimor_star *= pr;
    pry_dimor_star *= pr;
    prz_dimor_star *= pr;
    
    
#pragma omp parallel for
    for (int i=0;i<N;i++){
        
        struct reb_particle p = particles[i];             // cache
        if ( p.m > 0. ) continue;                         // Only dust particles feel radiation forces
        
        // set up vector
        const double prx_star  = p.x-star.x;
        const double pry_star  = p.y-star.y;
        const double prz_star  = p.z-star.z;
        const double prx_didy = p.x-Didymos.x;
        const double pry_didy = p.y-Didymos.y;
        const double prz_didy = p.z-Didymos.z;
        const double prx_dimor = p.x-Dimorphos.x;
        const double pry_dimor = p.y-Dimorphos.y;
        const double prz_dimor = p.z-Dimorphos.z;
        double dfactor;
        unsigned int flag_noshadow = 1;
        
        /* radiation force */
        // check if is in the shadow of Didymos
        pr = prx_didy*prx_didy_star + pry_didy*pry_didy_star + prz_didy*prz_didy_star;
        if ( pr > 0.0 ) {
            double prx_body_rad = prx_didy - pr*prx_didy_star;  // radial vector of dust relative to the Sun-Didymos direction
            double pry_body_rad = pry_didy - pr*pry_didy_star;  // radial vector of dust relative to the Sun-Didymos direction
            double prz_body_rad = prz_didy - pr*prz_didy_star;  // radial vector of dust relative to the Sun-Didymos direction
            if ( prx_body_rad*prx_body_rad + pry_body_rad*pry_body_rad + prz_body_rad*prz_body_rad < Rsq_didy )
                flag_noshadow = 0;
        }
        
        // check if is in the shadow of Dimorphos
        pr = prx_dimor*prx_dimor_star + pry_dimor*pry_dimor_star + prz_dimor*prz_dimor_star;
        if ( pr > 0.0 ) {
            double prx_body_rad = prx_dimor - pr*prx_dimor_star;  // radial vector of dust relative to the Sun-Dimorphos direction
            double pry_body_rad = pry_dimor - pr*pry_dimor_star;  // radial vector of dust relative to the Sun-Dimorphos direction
            double prz_body_rad = prz_dimor - pr*prz_dimor_star;  // radial vector of dust relative to the Sun-Dimorphos direction
            if ( prx_body_rad*prx_body_rad + pry_body_rad*pry_body_rad + prz_body_rad*prz_body_rad < Rsq_dimor )
                flag_noshadow = 0;
        }
        
        // add radiation pressure if not in shadow
        if ( flag_noshadow ) {
            pr = sqrt(prx_star*prx_star + pry_star*pry_star + prz_star*prz_star);     // distance relative to star
            //const double prvx = p.vx-star.vx;
            //const double prvy = p.vy-star.vy;
            //const double prvz = p.vz-star.vz;
            //const double rdot = (prvx*prx_star + prvy*pry_star + prvz*prz_star)/pr;     // radial velocity relative to star
            dfactor = SRP_coe/r_dust * pow(AU/pr,2.0);

            // Equation (5) of Burns, Lamy, Soter (1979)
            //particles[i].ax += dfactor*((1.-rdot/c)*prx_star/pr - prvx/c);
            //particles[i].ay += dfactor*((1.-rdot/c)*pry_star/pr - prvy/c);
            //particles[i].az += dfactor*((1.-rdot/c)*prz_star/pr - prvz/c);
            particles[i].ax += dfactor*prx_star/pr;
            particles[i].ay += dfactor*pry_star/pr;
            particles[i].az += dfactor*prz_star/pr;
        }
        
        // J2 of Didymos
        pr   = prx_didy*prx_didy + pry_didy*pry_didy + prz_didy*prz_didy;
        dfactor  = 3.0*r->G*J2_didy*Didymos.m*Didymos.r*Didymos.r/2./pow(pr,3.5);
        particles[i].ax += dfactor*prx_didy*(prx_didy*prx_didy + pry_didy*pry_didy - 4.*prz_didy*prz_didy);
        particles[i].ay += dfactor*pry_didy*(prx_didy*prx_didy + pry_didy*pry_didy - 4.*prz_didy*prz_didy);
        particles[i].az += dfactor*prz_didy*(3.*(prx_didy*prx_didy + pry_didy*pry_didy) - 2.*prz_didy*prz_didy);
        
        // J2 of Dimorphos
        pr   = prx_dimor*prx_dimor + pry_dimor*pry_dimor + prz_dimor*prz_dimor;
        dfactor  = 3.0*r->G*J2_dimor*Dimorphos.m*Dimorphos.r*Dimorphos.r/2./pow(pr,3.5);
        particles[i].ax += dfactor*prx_dimor*(prx_dimor*prx_dimor + pry_dimor*pry_dimor - 4.*prz_dimor*prz_dimor);
        particles[i].ay += dfactor*pry_dimor*(prx_dimor*prx_dimor + pry_dimor*pry_dimor - 4.*prz_dimor*prz_dimor);
        particles[i].az += dfactor*prz_dimor*(3.*(prx_dimor*prx_dimor + pry_dimor*pry_dimor) - 2.*prz_dimor*prz_dimor);
    }
}


void reb_move_to_Didymos(struct reb_simulation* const r){
    const int N_real = r->N - r->N_var;
    if (N_real>0){
        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle hel = r->particles[0];
        // Note: Variational particles will not be affected.
        for (int i=1;i<N_real;i++){
            particles[i].x  -= hel.x;
            particles[i].y  -= hel.y;
            particles[i].z  -= hel.z;
        }
        r->particles[0].x = 0.;
        r->particles[0].y = 0.;
        r->particles[0].z = 0.;
    }
}


void heartbeat(struct reb_simulation* r){
    
    // remove collide and escaped particles
    if(reb_simulation_output_check(r, 60.0)){
        
        struct reb_particle* particles = r->particles;
        const struct reb_particle Didymos = particles[0];
        const struct reb_particle Dimorphos = particles[1];
        int N = r->N;
        struct reb_vec3d vDis;
        double dDisSQ_Didy;
        unsigned int N_remove = 0;
        unsigned int flag_remove;
        
        // delete and record collided particles
        FILE* f_c = fopen("collide.txt","ab+");
        if ( f_c == NULL){
            reb_simulation_error(r, "Can not open file: collide.txt.");
            return;
        }

        for ( int i=0;i<N;i++ ) {
            
            const struct reb_particle p = particles[i-N_remove];       // cache
            if ( p.m > 0. ) continue;                         // Only delete dust particles
            
            flag_remove = 0;
            vDis.x = p.x-Didymos.x;
            vDis.y = p.y-Didymos.y;
            vDis.z = p.z-Didymos.z;
            dDisSQ_Didy = vDis.x*vDis.x + vDis.y*vDis.y + vDis.z*vDis.z;
            
            if ( dDisSQ_Didy < Rsq_didy )
                flag_remove = 1; // collide with Didymos
            else if ( ( (p.x-Dimorphos.x)*(p.x-Dimorphos.x) + (p.y-Dimorphos.y)*(p.y-Dimorphos.y) + (p.z-Dimorphos.z)*(p.z-Dimorphos.z) ) < Rsq_dimor )
                flag_remove = 2; // collide with Dimorphos
            else if ( dDisSQ_Didy > Rsq_hill )
                flag_remove = 3; // escape particles
                
            if ( flag_remove > 0 ) {
                
                fwrite( &(flag_remove), sizeof(int), 1, f_c );
                fwrite( &(p.hash), sizeof(int), 1, f_c );
                fwrite( &(r->t), sizeof(double), 1, f_c );
                
                reb_simulation_remove_particle( r, i-N_remove, 1 );
                N_remove++;
            }
        }
        fclose(f_c);
        
        //reb_move_to_hel(r);
        reb_move_to_Didymos(r);
    }
    
    //  output all particles
    if(reb_simulation_output_check(r, 3600.0)){
        
        struct reb_particle* particles = r->particles;
        const int N = r->N;
        double di;
        
        reb_simulation_output_timing(r, tmax);
        
        // output particle position and velocity
        FILE* fp = fopen("particles.txt","ab+");
        if ( fp == NULL){
            reb_simulation_error(r, "Can not open file: particles.txt.");
            return;
        }
        
        fwrite( &(N), sizeof(int), 1, fp);
        fwrite( &(r->t), sizeof(double), 1, fp);
        fwrite( &(r_dust), sizeof(double), 1, fp);
        for ( int i=0; i<N; i++ ) {
            const struct reb_particle p = particles[i];
            di = (double)p.hash;
            fwrite( &(di), sizeof(double), 1, fp);
            fwrite( &(p.x), sizeof(double), 1, fp);
            fwrite( &(p.y), sizeof(double), 1, fp);
            fwrite( &(p.z), sizeof(double), 1, fp);
            fwrite( &(p.vx), sizeof(double), 1, fp);
            fwrite( &(p.vy), sizeof(double), 1, fp);
            fwrite( &(p.vz), sizeof(double), 1, fp);
        }
        fclose(fp);
    }
}
