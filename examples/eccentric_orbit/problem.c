/**
 * Highly eccentric orbits
 *
 * This example uses the IAS15 integrator to simulate
 * a very eccentric planetary orbit. The integrator
 * automatically adjusts the timestep so that the pericentre passages
 * resolved with high accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double tmax;
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->G            = 1;        // Gravitational constant
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;

    double e_testparticle     = 1.-1e-7;    
    double mass_scale    = 1.;        // Some integrators have problems when changing the mass scale, IAS15 does not. 
    double size_scale    = 1;        // Some integrators have problems when changing the size scale, IAS15 does not.

    struct reb_particle star = {0}; 
    star.m  = mass_scale;
    reb_simulation_add(r, star); 
    
    struct reb_particle planet; 
    planet.m  = 0;
    planet.x  = size_scale*(1.-e_testparticle); 
    planet.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle)*mass_scale/size_scale);
    reb_simulation_add(r, planet); 
    
    reb_simulation_move_to_com(r);
    
    // initial timestep
    r->dt             = 1e-13*sqrt(size_scale*size_scale*size_scale/mass_scale); 
    tmax            = 1e2*2.*M_PI*sqrt(size_scale*size_scale*size_scale/mass_scale);

    reb_simulation_integrate(r, tmax);    
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r,tmax/10000.)){        // outputs to the screen
        reb_simulation_output_timing(r, tmax);
    }
    // Output the time and the current timestep. Plot it to see how IAS15 automatically reduces the timestep at pericentre. 
    FILE* of = fopen("timestep.txt","ab"); 
    fprintf(of,"%e\t%e\t\n",r->t/tmax,r->dt/tmax);
    fclose(of);
}

