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
    /* Several events can make particles move around in memory. This means the user should not assume they can always refer to 
     * the same particle through r->particles[index] or through a pointer to the particle that is set at the beginning of the 
     * simulation.  The reliable way to access particles is through hashes.*/

    struct reb_simulation* r = reb_simulation_create();

    struct reb_particle p = {0};
    p.m = 1.;
    p.hash = reb_hash("Sun");
    reb_simulation_add(r, p);

    printf("We can now reference the particle like this: m=%f, or like this: m=%f\n", r->particles[0].m, reb_simulation_particle_by_hash(r, reb_hash("Sun"))->m);

    for (int i=1; i <= 200; i++){
        struct reb_particle tp = {0};
        tp.x = i;       // put particles progressively farther away
        tp.hash = i;
        reb_simulation_add(r, tp);
        // ... initialize rest of particle variables.
    }

    /* The advantage of hashes is that if particles are ejected or otherwise removed from the simulation (or you are using 
     * the tree code), the indices in the particles array will get scrambled.  Accessing particles through their hash 
     * guarantees you get back the particle you intended.
     */

    printf("r->particles[200] hash=%u, x=%f\n", r->particles[200].hash, r->particles[200].x);

    reb_simulation_remove_particle_by_hash(r, reb_hash("Sun"), 0);

    /* The remove function has moved particles[200] to index 0:*/

    printf("r->particles[0] hash=%u, x=%f\n", r->particles[0].hash, r->particles[0].x);

    /* Rather than worry about the internals of what the remove function, the tree code etc. do, we can get it by hash. 
     * When you are not sure, assigning particles hashes and accessing them through them is always safe.
     * We can use reb_hash for a string, or just pass an unsigned integer we assigned directly.*/

    struct reb_particle* last = reb_simulation_particle_by_hash(r, 200);

    printf("Using hash: hash=%u, x=%f\n\n", last->hash, last->x); 

    /* Note that if the particle is not found in the simulation, reb_simulation_particle_by_hash returns a NULL pointer.
     * This allows you to check if particles are still in the simulation when you don't know ahead of time, but means that if you 
     * mistakenly access a removed particle, you'll get a segmentation fault.*/

    struct reb_particle* sunptr = reb_simulation_particle_by_hash(r, reb_hash("Sun"));

    if (sunptr == NULL){
        printf("Whoops!  Already removed particle.\n");
    }
    else{
        printf("Mass = %f\n\n", sunptr->m); // would cause segmentation fault in this case
    }

    /* The user is responsible for making sure the hashes don't clash. If two particles share the same hash, reb_simulation_particle_by_hash
     * could return either particle.  2 hashes generated with the reb_hash hash function have a ~1e-9 chance of clashing.
     * The most common case is assigning a hash of 0:
     */

    reb_simulation_free(r);
    r = reb_simulation_create();

    struct reb_particle sun = {0};
    sun.m = 1.;
    sun.hash = 0;
    reb_simulation_add(r, sun);

    struct reb_particle earth = {0};
    earth.x = 1.;
    earth.vy = 1.;
    reb_simulation_add(r, earth);

    printf("Sun's x position = %f\n", reb_simulation_particle_by_hash(r, 0)->x);

    /* The above line prints x=1 for the Sun's x position, which is not what we wanted.  The problem is we also set earth's hash to 0
     * when we initialized the structure to {0}!  We can use the 0 hash as long as we make sure we assign the hashes of all particles in the simulation:
     */

    r->particles[1].hash = reb_hash("earth");
    printf("Sun's x position = %f\n", reb_simulation_particle_by_hash(r, 0)->x);


    reb_simulation_free(r);
}

