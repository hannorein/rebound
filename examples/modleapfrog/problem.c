/**
 * Regularized Leap-frog
 *
 * This integrator is good at resolving highly eccentric orbits and 
 * close encounters of massive bodies. In this simulation, you
 * can visualize the increase in time resolution near periapsis as 
 * the simulation progresses through time nonlinearly.
 */
#include "rebound.h"
#include <math.h>


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
    r->usleep = 40000; //May need to adjust per device

    // Setup system
    reb_simulation_set_integrator(r, "modleapfrog");
    double M = 1000;    //primary mass
    double m = 1;       //secondary mass
    double a = 1;       //secondary semi-major axis
    double e = 0.99;    //secondary eccentricity
    double N = 100;     //target steps/orbit

    struct reb_particle primary = {0};
    primary.m = M;
    reb_simulation_add(r, primary);

    struct reb_particle secondary = reb_particle_from_orbit(r->G, primary, m, a, e, 0.0, 0.0, 0.0, 0.0);
    reb_simulation_add(r, secondary);
    reb_simulation_move_to_com(r);    

    r->dt =  2* M_PI * m * sqrt(M * a/ r->G)/ N; //Formula for fictional timestep dtau = 2pi*m*sqrt(Ma/G)/N

    reb_simulation_integrate(r, INFINITY);
}
