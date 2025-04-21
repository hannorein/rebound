#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct reb_simulation* setup(){
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_start_server(r, 1234);
    r->rand_seed = 0;
    reb_simulation_add_fmt(r, "m", 1.0);              
    //reb_simulation_add_fmt(r, "m", 0.32);              
    // reb_simulation_add_fmt(r, "m a e", 5e-5, 2.0, 0.01 );
    // reb_simulation_add_fmt(r, "m a e inc", 10e-3, 50.0, 0.52, 80.0/180.0*M_PI);
    for (int i=0; i<100;i++){
//        reb_simulation_add_fmt(r, "m a e E omega inc", 1e-6, 1.0+0.1*reb_random_normal(r,1.0), 0.94, reb_random_uniform(r,0,M_PI*2.0), reb_random_uniform(r,0,M_PI*2.0), 0.01*reb_random_normal(r, 1.0));
    }
    for (int i=0; i<10;i++){
        reb_simulation_add_fmt(r, "m r a e E inc", 1e-6, 0.1, 1.+0.1*reb_random_normal(r,1.0), 0.05, reb_random_uniform(r,0,M_PI*2.0), 0.01*reb_random_normal(r, 1.0));
    }
    reb_simulation_move_to_com(r);
 
    //struct reb_orbit o = reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);
    //r->dt = o.P/20.0;
    r->dt = 2.0*M_PI*1e-2;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    return r;
}

void test(struct reb_simulation* r, char* name){
    printf("%-20s:", name);
    fflush(stdout);
    double E0 = reb_simulation_energy(r);
    r->exact_finish_time = 0;
    //reb_simulation_integrate(r,2.0*M_PI*1e6);
    reb_simulation_integrate(r,100.);
    double E1 = reb_simulation_energy(r);

    printf(" dE = %e   runtime %5.1f s   N = %d\n", fabs((E0-E1)/E0), r->walltime, r->N);

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
    
    r = setup();
    r->integrator = REB_INTEGRATOR_BRACE;
    r->ri_brace.encounter_integrator = REB_BRACE_EC_IAS15;
    test(r, "brace-ias15");

}

