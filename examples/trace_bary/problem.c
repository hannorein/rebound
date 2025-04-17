#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct reb_simulation* setup(){
    struct reb_simulation* r = reb_simulation_create();
    
    //reb_simulation_start_server(r, 1234);

    reb_simulation_add_fmt(r, "m", 1.);              
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1);
    reb_simulation_add_fmt(r, "m a e", 1e-3, 2.4, 0.1);
    reb_simulation_move_to_com(r);
    return r;
}

void test(struct reb_simulation* r, char* name){
    double E0 = reb_simulation_energy(r);
    r->exact_finish_time = 0;
    reb_simulation_integrate(r,100.);
    double E1 = reb_simulation_energy(r);

    printf("%-20s: dE = %e \n", name, fabs((E0-E1)/E0));

    // Cleanup 
    reb_simulation_free(r);
}

int main(int argc, char* argv[]) {
    
    struct reb_simulation* r;
    
    r = setup();
    test(r, "ias15");
    
    
    r = setup();
    r->dt = 1e-1;
    r->integrator = REB_INTEGRATOR_WHFAST;
    test(r, "whfast");
    
    
    r = setup();
    r->dt = 1e-1;
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    test(r, "whfast-dh");
    
    r = setup();
    r->dt = 1e-1;
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_BARYCENTRIC;
    test(r, "whfast-bary");
    
    r = setup();
    r->dt = 1e-1;
    r->integrator = REB_INTEGRATOR_TRACE;
    test(r, "trace");

    r = setup();
    r->dt = 1e-1;
    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_trace.coordinates = REB_TRACE_COORDINATES_BARYCENTRIC;
    test(r, "trace-bary");

}

