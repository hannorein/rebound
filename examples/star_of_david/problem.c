/**
 * Star of David
 * 
 * This example uses the IAS15 integrator
 * to integrate the "Star od David", a four body system consisting of two
 * binaries orbiting each other. Note that the time is running backwards,
 * which illustrates that IAS15 can handle both forward and backward in time
 * integrations. The initial conditions are by Robert Vanderbei.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    r->integrator = REB_INTEGRATOR_IAS15;
    r->dt = -1;
    r->usleep = 10000;   // Slowing down integrator (for visualization only)

    struct reb_particle p = {0};
    p.m = 1.;
    p.z = 0.;
    p.vz = 0.;
    
    p.x =  -1.842389804706855; p.y =  -1.063801316823613; 
    p.vx =  -0.012073765486548; p.vy =   0.021537467220014; 
    reb_simulation_add(r, p);

    p.x =  -0.689515464218133; p.y =  -0.398759403276399; 
    p.vx =   0.637331229856386; p.vy =  -1.103822313621890; 
    reb_simulation_add(r, p);
    
    p.x =   0.689515464218133; p.y =   0.398759403276399; 
    p.vx =  -0.637331229856386; p.vy =   1.103822313621890; 
    reb_simulation_add(r, p);
    
    p.x =   1.842389804706855; p.y =   1.063801316823613; 
    p.vx =   0.012073765486548; p.vy =  -0.021537467220014; 
    reb_simulation_add(r, p);
    
    reb_simulation_move_to_com(r);

    reb_simulation_integrate(r, INFINITY);

    reb_simulation_free(r);
}
