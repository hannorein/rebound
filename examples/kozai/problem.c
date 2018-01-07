/**
 * Kozai cycles
 * 
 * This example uses the IAS15 integrator to simulate
 * a Lidov Kozai cycle of a planet perturbed by a distant star. 
 * The integrator automatically adjusts the timestep so that 
 * even very high eccentricity encounters are resolved with high 
 * accuracy.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double tmax = 1.6e4;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt             = M_PI*1e-2;     // initial timestep
    r->integrator        = REB_INTEGRATOR_IAS15;
    r->heartbeat        = heartbeat;

    // Initial conditions
    
    struct reb_particle star = {0}; 
    star.m  = 1;
    reb_add(r, star); 
    
    // The planet (a zero mass test particle)
    struct reb_particle planet = {0}; 
    double e_testparticle = 0;
    planet.m  = 0.;
    planet.x  = 1.-e_testparticle;
    planet.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle));
    reb_add(r, planet); 
    
    // The perturber
    struct reb_particle perturber = {0}; 
    perturber.x  = 10; 
    double inc_perturber = 89.9;
    perturber.m  = 1;
    perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt((star.m+perturber.m)/perturber.x); 
    perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt((star.m+perturber.m)/perturber.x); 
    reb_add(r, perturber); 

    reb_move_to_com(r);
    
    system("rm -v orbits.txt");        // delete previous output file

    reb_integrate(r, tmax);
}

void heartbeat(struct reb_simulation* r){
    if(reb_output_check(r, 20.*M_PI)){        // outputs to the screen
        reb_output_timing(r, tmax);
    }
    if(reb_output_check(r, 12.)){            // outputs to a file
        reb_output_orbits(r, "orbits.txt");
    }
}
