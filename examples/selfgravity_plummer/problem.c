/**
 * A self-gravitating Plummer sphere
 *
 * A self-gravitating Plummer sphere is integrated using
 * the leap frog integrator. Collisions are not resolved. Note that the
 * fixed timestep might not allow you to resolve individual two-body
 * encounters. An alternative integrator is IAS15 which
 * comes with adaptive timestepping.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup system characteristics
    int _N = 100;            // Number of particles
    double G = 1;            // Gravitational constant
    double M = 1;            // Total mass of the cluster
    double R = 1;            // Radius of the cluster
    double E = 3./64.*M_PI*M*M/R;   // Energy of the cluster
    double r0 = 16./(3.*M_PI)*R;    // Chacateristic length scale
    double t0 = r->G*pow(M,5./2.)*pow(4.*E,-3./2.)*(double)_N/log(0.4*(double)_N); // Rellaxation time
    printf("Characteristic size:              %f\n", r0);
    printf("Characteristic time (relaxation): %f\n", t0);

    // Setup constants
    r->G         = G;        
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->dt         = 2e-5*t0;     // timestep
    r->softening     = 0.01*r0;    // Softening parameter
    r->heartbeat    = heartbeat;
    
    reb_simulation_configure_box(r, 20.*r0, 1, 1, 1);
    reb_simulation_add_plummer(r, _N, M, R);    // Adds particles
    reb_simulation_move_to_com(r); 
    reb_simulation_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.0*r->dt)){
        reb_simulation_output_timing(r, 0);
    }
}
