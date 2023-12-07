/**
 * Velocity dependent drag force
 *
 * This is a very simple example on how to implement a velocity 
 * dependent drag force. The example uses the IAS15 integrator, which 
 * is ideally suited to handle non-conservative forces.
 * No gravitational forces or collisions are present.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);

double tmax = 40.;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
   
    // Setup constants
    r->dt             = 1e-4;        // initial timestep.
    r->integrator     = REB_INTEGRATOR_IAS15;
    r->gravity        = REB_GRAVITY_NONE;

    // Setup callback function for velocity dependent forces.
    r->additional_forces = additional_forces;
    r->force_is_velocity_dependent = 1;
    // Setup callback function for outputs.
    r->heartbeat   = heartbeat;
    // Slow down integration (for visualization only)
    r->usleep = 10000;  
    
    struct reb_particle p = {0}; 
    p.m      = 0;    // massless
    p.x     = 1;
    p.vx     = -1;
    reb_simulation_add(r, p); 

    // Delete previous output
    remove("r.txt");    

    // Do the integration
    reb_simulation_integrate(r, tmax);
}

void additional_forces(struct reb_simulation* const r){
    // Simplest velocity dependent drag force.
    double dragcoefficient = 1;
    struct reb_particle* const particles = r->particles;
    const int N = r->N;
    for (int i=0;i<N;i++){
        particles[i].ax = -dragcoefficient*particles[i].vx;
        particles[i].ay = -dragcoefficient*particles[i].vy;
        particles[i].az = -dragcoefficient*particles[i].vz;
    }
}

void heartbeat(struct reb_simulation* const r){
    // Output some information to the screen every 100th timestep
    if(reb_simulation_output_check(r, 100.*r->dt)){
        reb_simulation_output_timing(r, tmax);
    }
    // Output the particle position to a file every timestep.
    const struct reb_particle* const particles = r->particles;
    FILE* f = fopen("r.txt","ab");
    fprintf(f,"%e\t%e\t%e\n",r->t,particles[0].x, particles[1].vx);
    fclose(f);
}
