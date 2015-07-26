/**
 * How to use unique ids to identify particles
 *
 * This example shows how to assign ids to particles, and demonstrates different 
 * options for removing particles from the simulation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void print_ids(struct reb_simulation* r){
	printf("ids = ");
	for (int i=0;i<r->N;i++){
		printf("%d ", r->particles[i].id);
	}
	printf("\n");
}


int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Just demonstrating ids, so set initial conditions arbitrarily (and ids in the order they are added)
	// ids can be set to any integer
	for (int i=0;i<10;i++){
		struct reb_particle p = {0};
		p.id = i; 
		reb_add(r, p); 
	}

	printf("Initial ids:\n");
	print_ids(r);

	int success;
	int keepSorted = 0;
	printf("\nTry to remove index 3...\n");
	success = reb_remove(r, 3, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_ids(r);
	printf("Because keepSorted = 0, last particle replaced removed particle and indices got scrambled:\n\n");

	keepSorted = 1;
	printf("Try to remove index 6 while preserving the order with keepSorted=1...\n");
	success = reb_remove(r, 6, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_ids(r);

	printf("\nWe can also remove particles by id.  Try to remove id=5...\n");
	success = reb_remove_by_id(r, 5, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	print_ids(r);
	
	printf("\nIf we try to remove an index > N or an id that doesn't exist, we get a warning and no particle is removed:\n");
	printf("Try to remove index 15...\n");
	success = reb_remove(r, 15, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
	printf("Try to remove id=3...\n");
	success = reb_remove_by_id(r, 3, keepSorted);
	if (success){
		printf("Particle successfully removed\n");
	}
}

