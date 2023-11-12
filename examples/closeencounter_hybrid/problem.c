/**
 * Close encounters with MERCURIUS integrator
 * 
 * This example integrates a densely packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. 
 * This is a test case for the MERCURIUS integrator.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    r->dt = 0.0012*2.*M_PI;                
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->ri_mercurius.r_crit_hill = 3;            // By default the switching radius is three times the hill radius
    r->heartbeat = heartbeat;

    struct reb_particle star = {0};
    star.m = 1;
    reb_simulation_add(r, star);
    
    // Add planets
    int N_planets = 3;
    for (int i=0;i<N_planets;i++){
        double a = 1.+.1*(double)i;        // semi major axis
        double v = sqrt(1./a);             // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 2e-3; 
        planet.x = a; 
        planet.vy = v;
        reb_simulation_add(r, planet); 
    }
    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_simulation_energy(r);
    remove("energy.txt");

    reb_simulation_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.*2.*M_PI)){  
        reb_simulation_output_timing(r, 0);
    }
    if (reb_simulation_output_check(r, 2.*M_PI)){  
        // Once per year, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","ab");
        reb_simulation_synchronize(r);
        double e = reb_simulation_energy(r);
        fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        fclose(f);
    }
}

