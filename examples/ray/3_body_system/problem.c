#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void add_central_particle(struct reb_simulation* const r){
    struct reb_particle p = {0};
    p.m = 1.0;
    reb_simulation_add(r, p);
}

static void add_circular_particle(struct reb_simulation* const r, const double mass, const double radius, const double f){
    struct reb_particle p = {0};
    p.m = mass;
    p.x = radius*cos(f);
    p.y = radius*sin(f);
    const double speed = sqrt(1.0/radius);
    p.vx = -speed*sin(f);
    p.vy = speed*cos(f);
    reb_simulation_add(r, p);
}

static void add_three_body_system(struct reb_simulation* const r){
    r->G = 1.0;

    add_central_particle(r);

    const double m = 5.e-8;
    add_circular_particle(r, m, 1.0, 0.0);
    add_circular_particle(r, m, 2.0, 0.0);
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    r->dt = 0.1;
    reb_simulation_set_integrator(r, "whfast_hj");
    // reb_simulation_set_integrator(r, "whfast");  

    add_three_body_system(r);
    reb_simulation_move_to_com(r);

    const double E0 = reb_simulation_energy(r);
    reb_simulation_integrate(r, 100.*M_PI*2.);
    const double E1 = reb_simulation_energy(r);

    printf("dE = %e\n", fabs((E0 - E1)/E0));

    reb_simulation_free(r);
}
// dE = 1.982657e-11 whfast
// dE = 1.398063e-11 whfast_hj