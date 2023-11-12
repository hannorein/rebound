/**
 * Radiation forces
 *
 * This example provides an implementation of the 
 * Poynting-Robertson effect. The code is using the IAS15 integrator
 * which is ideally suited for this velocity dependent force.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void radiation_forces(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

double betaparticles = 0.01;     // beta parameter, defined as the ratio of radiation pressure over gravity
double tmax = 1e5;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // setup constants
    r->dt                           = 1e-3;             // initial timestep
    r->integrator                   = REB_INTEGRATOR_IAS15;
    r->ri_ias15.epsilon             = 1e-4;             // accuracy parameter
    r->N_active                     = 1;                // the star is the only massive particle
    r->force_is_velocity_dependent  = 1;
    r->additional_forces            = radiation_forces; // setup callback function for velocity dependent forces
    r->heartbeat                    = heartbeat;
    
    // star is at rest at origin
    struct reb_particle star = {0};
    star.m  = 1.;
    reb_simulation_add(r, star);

    // dust particles are initially on a circular orbit
    while(r->N<2){
        struct reb_particle p = {0}; 
        p.m  = 0;                    // massless
        double a = 1.;                    // a = 1 AU
        double v = sqrt(r->G*(star.m*(1.-betaparticles))/a);
        double phi = reb_random_uniform(r, 0,2.*M_PI);        // random phase
        p.x  = a*sin(phi);  p.y  = a*cos(phi); 
        p.vx = -v*cos(phi); p.vy = v*sin(phi);
        reb_simulation_add(r, p); 
    }

    remove("radius.txt");                    // remove previous output

    reb_simulation_integrate(r, tmax);
}

void radiation_forces(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    const int N = r->N;
    const struct reb_particle star = particles[0];    // cache
#pragma omp parallel for
    for (int i=0;i<N;i++){
        const struct reb_particle p = particles[i];   // cache
        if (p.m!=0.) continue;                        // only dust particles feel radiation forces
        const double prx  = p.x-star.x;
        const double pry  = p.y-star.y;
        const double prz  = p.z-star.z;
        const double pr   = sqrt(prx*prx + pry*pry + prz*prz);         // distance relative to star
        const double prvx = p.vx-star.vx;
        const double prvy = p.vy-star.vy;
        const double prvz = p.vz-star.vz;

        const double c         = 1.006491504759635e+04;         // speed of light in unit of G=1, M_sun=1, 1year=1
        const double rdot     = (prvx*prx + prvy*pry + prvz*prz)/pr;     // radial velocity relative to star
        const double F_r     = betaparticles*r->G*star.m/(pr*pr);

        // Equation (5) of Burns, Lamy, Soter (1979)
        particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
        particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
        particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);
    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r, 400.)){                        // print some information to screen
        reb_simulation_output_timing(r, tmax);;
    }
    if(reb_simulation_output_check(r, M_PI*2.*1000.)){                     // output radial distance every 1000 years
        FILE* f = fopen("radius.txt","ab");
        struct reb_particle* particles = r->particles;
        const struct reb_particle star = particles[0];
        const int N = r->N;
        for (int i=1;i<N;i++){
            const struct reb_particle p = particles[i]; 
            const double prx  = p.x-star.x;
            const double pry  = p.y-star.y;
            const double prz  = p.z-star.z;
            const double pr   = sqrt(prx*prx + pry*pry + prz*prz);     // distance relative to star
            fprintf(f,"%e\t%e\n",r->t,pr);
        }
        fclose(f);
    }
}
