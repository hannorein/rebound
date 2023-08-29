/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m", 1.);                // Central object
    reb_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle 

    reb_integrate(r,100.);

    struct reb_simulation* r2 = reb_create_simulation();

    reb_output_binary(r2, "test");

    struct reb_simulation* r3 = reb_create_simulation_from_binary("test");
    reb_diff_simulations(r3, r,1);

//    reb_free_simulation(r);
  //  reb_free_simulation(r2);
}

