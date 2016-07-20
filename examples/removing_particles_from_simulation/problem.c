/**
 * Removing particles from simulations.
 *
 * This example demonstrates different 
 * options for removing particles from the simulation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void print_hashes(struct reb_simulation* r){
	printf("hashes = ");
	for (int i=0;i<r->N;i++){
		printf("%u ", r->particles[i].hash);
	}
	printf("\n");
}


int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();

	for (int i=0;i<9;i++){
        struct reb_particle p = {0};
        p.hash = i;
        reb_add(r, p);
	}

    struct reb_particle p = {0};
    p.hash = reb_hash("Planet 9");
    reb_add(r, p);

    printf("Initial hashes:\n");
	print_hashes(r);

	int success;
	int keepSorted = 0;
	printf("\nTry to remove index 3 (4th particle)...\n");
	success = reb_remove(r, 3, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_hashes(r);
	printf("Because keepSorted = 0, last particle replaced removed particle and indices got scrambled:\n\n");

	keepSorted = 1;
	printf("Try to remove index 6 (7th particle)  while preserving the order with keepSorted=1...\n");
	success = reb_remove(r, 6, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_hashes(r);
	
    printf("\nWe can also remove particles by the hashes we assign them (this is robust to particles switching indices in the particles array during the simulation).\n");  
    printf("Try to remove Planet 9...\n");
    success = reb_remove_by_hash(r, reb_hash("Planet 9"), keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_hashes(r);
   
    printf("\nFinally, we can remove particles by their hash directly.\n");
	success = reb_remove_by_hash(r, 1, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_hashes(r);
	
	printf("\nAlso, if we try to remove an index > N, we get an error and no particle is removed:\n");
	printf("Try to remove index 15...\n");
	success = reb_remove(r, 15, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
}

