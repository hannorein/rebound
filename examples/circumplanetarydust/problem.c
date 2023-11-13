/**
 * Radiation forces on circumplanetary dust
 * 
 * This example shows how to integrate circumplanetary
 * dust particles using the IAS15 integrator.
 * The example sets the function pointer `additional_forces`
 * to a function that describes the radiation forces.
 * The example uses a beta parameter of 0.01. 
 * The output is custom too, outputting the semi-major axis of 
 * every dust particle relative to the planet. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void force_radiation(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

double betaparticles = 0.01;     // Beta parameter, defined as the ratio of radiation pressure over gravity.
double tmax = 1e6;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    // Setup constants
    r->integrator           = REB_INTEGRATOR_IAS15;
    r->dt                   = 1e-4;      // Initial timestep.
    r->N_active             = 2;         // Only the star and the planet are massive.
    r->additional_forces    = force_radiation;
    r->heartbeat            = heartbeat;
    r->usleep               = 5000;      // Slow down integration (for visualization only)
    
    // Star
    struct reb_particle star = {0};
    star.m  = 1.;
    reb_simulation_add(r, star);


    // planet 
    struct reb_particle planet = {0};
    planet.m  = 1e-3;
    planet.x  = 1; 
    planet.vy = sqrt(r->G*(star.m+planet.m)/planet.x);
    reb_simulation_add(r, planet);
    
    

    // Dust particles
    while(r->N<3){          // Three particles in total (star, planet, dust particle) 
        struct reb_particle p = {0}; 
        p.m  = 0;           // massless
        double _r = 0.001;  // distance from planet planet
        double v = sqrt(r->G*planet.m/_r);
        p.x  = _r;
        p.vy = v;
        p.x += planet.x;     p.y += planet.y;     p.z += planet.z;
        p.vx += planet.vx;     p.vy += planet.vy;     p.vz += planet.vz;
        reb_simulation_add(r, p); 
    }
    
    reb_simulation_move_to_com(r);

    remove("a.txt");    

    reb_simulation_integrate(r, tmax);
}

void force_radiation(struct reb_simulation* r){
    struct reb_particle* particles = r->particles;
    const struct reb_particle star = particles[0];            // cache
    const int N = r->N;
    const double G = r->G;
#pragma omp parallel for
    for (int i=0;i<N;i++){
        const struct reb_particle p = particles[i];             // cache
        if (p.m!=0.) continue;                         // Only dust particles feel radiation forces
        const double prx  = p.x-star.x;
        const double pry  = p.y-star.y;
        const double prz  = p.z-star.z;
        const double pr   = sqrt(prx*prx + pry*pry + prz*prz);     // distance relative to star
        const double prvx = p.vx-star.vx;
        const double prvy = p.vy-star.vy;
        const double prvz = p.vz-star.vz;

        const double c         = 1.006491504759635e+04;         // speed of light.
        const double rdot     = (prvx*prx + prvy*pry + prvz*prz)/pr;     // radial velocity relative to star
        const double F_r     = betaparticles*G*star.m/(pr*pr);

        // Equation (5) of Burns, Lamy, Soter (1979)
        particles[i].ax += F_r*((1.-rdot/c)*prx/pr - prvx/c);
        particles[i].ay += F_r*((1.-rdot/c)*pry/pr - prvy/c);
        particles[i].az += F_r*((1.-rdot/c)*prz/pr - prvz/c);
    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r, M_PI*2.)){
        reb_simulation_output_timing(r, tmax);
    }
    if(reb_simulation_output_check(r, M_PI*2.)){ // output every year
        FILE* f = fopen("a.txt","ab");
        const struct reb_particle* particles = r->particles;
        const struct reb_particle planet = particles[1];
        const double G = r->G;
        const double t = r->t;
        const int N = r->N;
        for (int i=2;i<N;i++){
            const struct reb_particle p = particles[i]; 
            const double prx  = p.x-planet.x;
            const double pry  = p.y-planet.y;
            const double prz  = p.z-planet.z;
            const double pr   = sqrt(prx*prx + pry*pry + prz*prz);     // distance relative to star
            
            const double pvx  = p.vx-planet.vx;
            const double pvy  = p.vy-planet.vy;
            const double pvz  = p.vz-planet.vz;
            const double pv   = sqrt(pvx*pvx + pvy*pvy + pvz*pvz);     // distance relative to star
            
            const double a = -G*planet.m/( pv*pv - 2.*G*planet.m/pr );            // semi major axis
            
            fprintf(f,"%e\t%e\n",t,a);
        }
        fclose(f);
    }
}
