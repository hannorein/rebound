#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_add_fmt(r, "m", 1.);
    reb_simulation_add_fmt(r, "m P", 1e-6, 1.0);
    reb_simulation_add_fmt(r, "m P", 1e-6, 2.0);
    reb_simulation_add_fmt(r, "m P", 1e-6, 3.0);
    reb_simulation_add_fmt(r, "m P", 1e-6, 4.0);
    reb_simulation_add_fmt(r, "m P", 1e-6, 5.0);
    reb_simulation_add_fmt(r, "m a e", 1e-6, -5.0, 2.0);
    reb_simulation_add_fmt(r, "m a e", 1e-6, -10.0, 2.0);
 //   reb_simulation_add_fmt(r, "m P", 10.0, 105.0);
//    reb_simulation_add_fmt(r, "m P primary", 1e-6, 1.0, r->particles[5]);
  //  reb_simulation_add_fmt(r, "m P primary", 1e-6, 0.01, r->particles[2]);
    
//    reb_simulation_add_fmt(r, "m", 1.);
//    reb_simulation_add_fmt(r, "m P primary", 0.01, 2.0, r->particles[0]);
//    reb_simulation_add_fmt(r, "m P primary", 0.01, 1.0, r->particles[0]);
//    reb_simulation_add_fmt(r, "m a e", 0.01, -1.0, 9.3);

//    reb_simulation_add_fmt(r, "m", 1.);
//    reb_simulation_add_fmt(r, "m P e", 1e-3, 9., 0.1);
//    reb_simulation_add_fmt(r, "m P e", 1e-3, 9., 0.1);
//    reb_simulation_add_fmt(r, "P e", 10.0, 0.1);
//    reb_simulation_add_fmt(r, "m P e primary", 1e-6, 4.0, 0.1, r->particles[1]);

    struct reb_orbit_hierarchy* oh = reb_orbit_hierarchy_create_from_simulation(r);
    reb_orbit_hierarchy_print(oh,r,0);
    printf("Orbit node is Jacobi: %d\n", reb_orbit_hierarchy_is_jacobi(oh));
    printf("Orbit node is Jacobi ordered: %d\n", reb_orbit_hierarchy_is_jacobi_ordered(oh, r));
    reb_orbit_hierarchy_free(oh);



    reb_simulation_free(r);
}

