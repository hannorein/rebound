#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt                   = 0.1;      
    reb_simulation_set_integrator(r, "whfast_hj");
    // reb_simulation_set_integrator(r, "whfast");

    // Initial conditions
    reb_simulation_add_fmt(r, "solar system"); // Built in dataset for testing.
    reb_simulation_move_to_com(r);

    double E0 = reb_simulation_energy(r);
    reb_simulation_integrate(r, 100*M_PI*2);
    double E1 = reb_simulation_energy(r);
    printf("dE = %e\n", fabs((E0-E1)/E0));
}

// dE = 4.464914e-10 whfast_hj before implementing
// dE = 4.464980e-10 whfast
// dE = 4.507626e-10 whfast_hj after implementing