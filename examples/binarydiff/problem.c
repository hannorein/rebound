/**
 * Comparing two binary files
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    r->integrator=REB_INTEGRATOR_IAS15;
    r->dt = 0.1;
    struct reb_particle p1 = {0}; 
    p1.m = 1.;
    reb_add(r, p1);  
    
    struct reb_particle p2 = {0};
    p2.x = 1;
    p2.vy = 1;
    reb_add(r, p2); 
    
    
    r->simulationarchive_filename = "sa.bin";
    r->simulationarchive_interval = 5;
    r->simulationarchive_version = 2;

    reb_integrate(r,100.);

    reb_output_binary(r, "s1.bin");
    reb_integrator_reset(r);
    r->integrator=REB_INTEGRATOR_WHFAST;
    struct reb_particle p3 = {0};
    p3.x = 2;
    p3.vy = 0.21;
    reb_add(r, p3); 
    reb_integrate(r,110.);
    reb_output_binary(r, "s2.bin");

    FILE* f1 = fopen("s1.bin","r");
    FILE* f2 = fopen("s2.bin","r");
    char* buf;
    size_t size;
    reb_binary_diff(f1,f2,&buf,&size);

}
