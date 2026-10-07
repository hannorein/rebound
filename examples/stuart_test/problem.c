
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the visualization web server.
    // Point your browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
    
    // Setup constants
    reb_simulation_set_integrator(r, "modleapfrog");

    r->usleep    = 1000;            // Slow down integration (for visualization only)
    r->dt = 1e-2;
    
    // Initial conditions
    struct reb_particle primary = {0};
    primary.m = 10;
    reb_simulation_add(r, primary);

    struct reb_particle secondary = reb_particle_from_orbit(r->G, primary, 0.1, 10, 0.999, 0.0, 0.0, 0.0, 2);
    reb_simulation_add(r, secondary);

    reb_simulation_move_to_com(r);

    r->exact_finish_time = 0;
    reb_simulation_integrate(r, 10.);

}

