/**
 * How to use hashes to identify particles
 *
 * This example shows how to assign hashes to particles
 * and how to access particles using hashes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    struct reb_particle p = {0};
    p.m = 1.;
    reb_simulation_add(r, p);
    reb_particle_set_name(&r->particles[0], "Sun");
    printf("registered names: %d\n", r->N_name_list);
    
    struct reb_particle p2 = {0};
    p2.m = 1.e-3;
    p2.name = "Earth";
    reb_simulation_add(r, p2);
    printf("registered names: %d\n", r->N_name_list);

    struct reb_particle* p3 = reb_simulation_get_particle_by_name(r, "Sun");
    if (p3){
        printf("Found p3. m=%f\n", p3->m);
    }else{
        printf("Did not find p3\n");
    }
    
    struct reb_particle* p4 = reb_simulation_get_particle_by_name(r, "Earth");
    if (p4){
        printf("Found p4. m=%f\n", p4->m);
    }else{
        printf("Did not find p4\n");
    }

    system("rm -rf out.bin");
    reb_simulation_save_to_file(r, "out.bin");
    reb_simulation_free(r);
    
    // Works by accident so far (memory of name_list never freed)
    r = reb_simulation_create_from_file("out.bin",-1);
    printf("registered names: %d\n", r->N_name_list);
    p3 = reb_simulation_get_particle_by_name(r, "Sun");
    printf("registered names: %d\n", r->N_name_list);
    if (p3){
        printf("Found p3. m=%f\n", p3->m);
    }else{
        printf("Did not find p3\n");
    }
    
    p4 = reb_simulation_get_particle_by_name(r, "Earth");
    if (p4){
        printf("Found p4. m=%f\n", p4->m);
    }else{
        printf("Did not find p4\n");
    }
    reb_simulation_free(r);
}

