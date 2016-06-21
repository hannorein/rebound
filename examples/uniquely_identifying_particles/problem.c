/**
 * How to use hashes to identify particles
 *
 * This example shows how to assign hashes to particles
 * and how to access particles using hashes.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void print_hashes(struct reb_simulation* r){
	printf("hashes = ");
	for (int i=0;i<r->N;i++){
		printf("%u ", r->particles[i].hash); // hashes are 32-bit unsigned integers
	}
	printf("\n");
}


int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
    struct reb_particle p = {0};
    p.m = 1.;
    reb_add(r, p);

    r->particles[0].hash = reb_tools_hash("Sun"); // We add a hash to particles[0] corresponding to the name "Sun"

    printf("We can now reference the particle like this: m=%f, or like this: m=%f\n", r->particles[0].m, reb_get_particle_by_name(r, "Sun")->m);
	
    /* The advantage of accessing particles with a unique hash over using the index is that the indices in the particles array can get scrambled
     * as particles are added/removed, if you're using the gravity tree code etc. Note reb_get_particle_by_name (and reb_get_particle_by_hash) below
     * return a pointer to the particle.  If you need to assign many hashes and don't want to give them
     * names, you can also have the simulation assign particles a unique hash: 
     */

    for (int i=1;i<10;i++){
		struct reb_particle p = {0};
		reb_add(r, p); 
        r->particles[i].hash = reb_generate_unique_hash(r); // this returns consecutive hashes, starting from a random number.
	}

	print_hashes(r);
    
    printf("We can also reference particles with their numeric hash directly.\n");
    uint32_t hash = reb_tools_hash("Sun"); 
    printf("Hash corresponding to the name Sun = %u\n", hash);
    printf("m of Sun = %f\n", reb_get_particle_by_hash(r,hash)->m);
    // You can therefore assign your own hash values manually, but if two particles share the same hash, get_particle_by_hash will return the first hit in particles array.
}

