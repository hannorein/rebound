#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_add_fmt(r, "m", 1.);                // Central object
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_simulation_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle 

    reb_simulation_set_integrator(r, reb_integrator_whfast);
    reb_simulation_integrate(r,1.);

    printf("r->particles   %zu\n", (r->particles));

    struct reb_simulation* r2 = reb_simulation_copy(r);
    printf("r->particles   %zu\n", (r->particles));
    printf("r2->particles  %zu\n", (r2->particles));

    reb_simulation_integrate(r,2.);
    reb_simulation_integrate(r2,2.);

    printf("r->particles   %zu\n", (r->particles));
    printf("r2->particles  %zu\n", (r2->particles));

    reb_simulation_free(r);
    reb_simulation_free(r2);
}

