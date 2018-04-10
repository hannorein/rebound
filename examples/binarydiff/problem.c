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

    reb_integrate(r,100.);

    reb_output_binary(r, "s1.bin");
    reb_integrator_reset(r);
    r->integrator=REB_INTEGRATOR_WHFAST;
    reb_integrate(r,110.);
    reb_output_binary(r, "s2.bin");

    FILE* f1 = fopen("s1.bin","r");
    FILE* f2 = fopen("s2.bin","r");
    FILE* diff = reb_binary_diff(f1,f2);
    fclose(diff);

}

