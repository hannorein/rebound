#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct reb_simulation* setup(){
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_start_server(r, 1234);
    r->rand_seed = 0;
    reb_simulation_add_fmt(r, "m", 1.);              
    for (int i=0; i<100;i++){
        reb_simulation_add_fmt(r, "m a e f", 1e-6, 1.+0.1*reb_random_normal(r,1.0), 0.1, reb_random_uniform(r,0,M_PI*2.0));
    }
    reb_simulation_move_to_com(r);
 
    r->dt = M_PI*2.0*1e-2;
    return r;
}

void test(struct reb_simulation* r, char* name){
    double E0 = reb_simulation_energy(r);
    r->exact_finish_time = 0;
    reb_simulation_integrate(r,100.);
    double E1 = reb_simulation_energy(r);

    printf("%-20s: dE = %e   runtime %.1f s\n", name, fabs((E0-E1)/E0), r->walltime);

    // Cleanup 
    reb_simulation_free(r);
}

int main(int argc, char* argv[]) {
    
    struct reb_simulation* r;
    
    r = setup();
    test(r, "ias15");
    
    r = setup();
    r->integrator = REB_INTEGRATOR_WHFAST;
    test(r, "whfast");
    
    r = setup();
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC;
    test(r, "whfast-dh");
    
    r = setup();
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_BARYCENTRIC;
    test(r, "whfast-bary");
    
    r = setup();
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    test(r, "mercurius");
    
    r = setup();
    r->integrator = REB_INTEGRATOR_TRACE;
    test(r, "trace");

    r = setup();
    r->integrator = REB_INTEGRATOR_BRACE;
    test(r, "brace");

}

