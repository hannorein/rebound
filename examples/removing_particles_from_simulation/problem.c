/**
 * Removing particles from simulations
 *
 * This example demonstrates different 
 * options for removing particles from the simulation.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_names(struct reb_simulation* r){
    for (int i=0;i<r->N;i++){
        printf("r->particles[%i].name = %s\n", i, r->particles[i].name);
    }
    printf("\n");
}


int main(int argc, char* argv[]){
    // Note that when using a tree (for gravity calculation or collision search), you need
    // to call reb_tree_update(r) after removing particles. Only then is the particle
    // removed. This helps avoiding to rebuild the tree multiple times if more than one
    // particle is removed at the same time.
    struct reb_simulation* r = reb_simulation_create();
    
    // The following lines add particles with names in three different ways.
    reb_simulation_add(r, (struct reb_particle){.name="Sun"});

    for (int i=1;i<9;i++){
        struct reb_particle p = {0};
        char name[256];
        sprintf(name, "Planet %d", i);
        p.name = name;
        reb_simulation_add(r, p);
    }

    struct reb_particle p = {.name = "Planet 9"};
    reb_simulation_add(r, p);

    printf("Initial names:\n");
    print_names(r);

    int error;
    int keep_sorted = 0;
    printf("\nTry to remove index 3 (Planet 2)...\n");
    error = reb_simulation_remove_particle(r, 3, keep_sorted);
    if (!error){
        printf("Particle successfully removed\n");
    }
    print_names(r);
    printf("Because keep_sorted = 0, last particle replaced the removed particle and indices got scrambled:\n\n");

    keep_sorted = 1;
    printf("Try to remove index 6 (Planet 7)  while preserving the order with keep_sorted=1...\n");
    error = reb_simulation_remove_particle(r, 6, keep_sorted);
    if (!error){
        printf("Particle successfully removed.\n");
    }
    print_names(r);
    
    printf("\nWe can also remove particles by the names we assign them (this is robust to particles switching indices in the particles array during the simulation).\n");  
    printf("Try to remove Planet 9...\n");
    error = reb_simulation_remove_particle_by_name(r, "Planet 9", keep_sorted);
    if (!error){
        printf("Particle successfully removed.\n");
    }
    print_names(r);
   
    printf("\nAlso, if we try to remove an index > N, we get an error and no particle is removed:\n");
    printf("Try to remove index 15...\n");
    error = reb_simulation_remove_particle(r, 15, keep_sorted);
    if (!error){
        printf("Particle successfully removed.\n");
    }else{
        printf("Particle not removed.\n");
    }

    reb_simulation_free(r);
}

