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

int main(int argc, char* argv[]){
    /* Several events can make particles move around in memory. This means the user should not assume they can always refer to 
     * the same particle through r->particles[index] or through a pointer to the particle that is set at the beginning of the 
     * simulation.  The reliable way to access particles is through hashes.*/

    struct reb_simulation* r = reb_create_simulation();

    /* In order for hashes to work correctly and faster, we have to add particles using reb_add_particle (rather than reb_add).
     * reb_add_particle adds a new particle to the simulation and returns a pointer to it.  It sets the particle simulation pointer
     * and hash to the default value UINT32_MAX but the particle is otherwise uninitialized. reb_add is kept for backwards compatibility.
     */

    struct reb_particle* p = reb_add_particle(r);
    p->m = 1.;
    p->hash = reb_tools_hash("Sun"); // We add a hash to particles[0] corresponding to the name "Sun"

    printf("We can now reference the particle like this: m=%f, or like this: m=%f\n", r->particles[0].m, reb_get_particle_by_hash(r, reb_tools_hash("Sun"))->m);

    for (int i=1; i <= 200; i++){
        struct reb_particle* tp = reb_add_particle(r);
        tp->hash = i;
        tp->x = i;      // put particles progressively farther away
        // ... initialize rest of particle variables.
    }

    /* One pitfall is that the particles array is reallocated in blocks of 128 particles.  If we try to update the Sun's mass with
     * our original pointer, it doesn't work because it's pointing to an old location in memory. */

    p->m = 1000.;
    printf("Sun mass still = %f after setting p->m = 1000.\n", r->particles[0].m);

    reb_get_particle_by_hash(r, reb_tools_hash("Sun"))->m = 1000.;
    printf("Sun mass correctly = %f using hash.\n\n", r->particles[0].m);

    /* If you always initialize the particle just after getting the pointer back from reb_add_particle there is no problem.*/
    
    /* The advantage of hashes is that if particles are ejected or otherwise removed from the simulation (or you are using 
     * the tree code), the indices in the particles array will get scrambled.  Accessing particles through their hash 
     * guarantees you get back the particle you intended.
     */

    printf("r->particles[200] hash=%u, x=%f\n", r->particles[200].hash, r->particles[200].x);

    reb_remove_by_hash(r, reb_tools_hash("Sun"), 0);

    /* The remove function has moved particles[200] to index 0:*/

    printf("r->particles[0] hash=%u, x=%f\n", r->particles[0].hash, r->particles[0].x);

    /* Rather than worry about the internals of what the remove function, the tree code etc. do, we can get it by hash. 
     * When you are not sure, assigning particles hashes and accessing them through them is always safe.
     * We can use reb_tools_hash for a string, or just pass an unsigned integer we assigned directly.*/

    struct reb_particle* last = reb_get_particle_by_hash(r, 200);

    printf("Using hash: hash=%u, x=%f\n\n", last->hash, last->x); 

    /* Note that if the particle is not found in the simulation, reb_get_particle_by_hash returns a NULL pointer.
     * This allows you to check if particles are still in the simulation when you don't know ahead of time, but means that if you 
     * mistakenly access a removed particle, you'll get a segmentation fault.*/

    struct reb_particle* sun = reb_get_particle_by_hash(r, reb_tools_hash("Sun"));

    if (sun == NULL){
        printf("Whoops!  Already removed particle.\n");
    }
    else{
        printf("Mass = %f\n", sun->m); // would cause segmentation fault in this case
    }

    /* The user is responsible for making sure the hashes don't clash. If two particles share the same hash, reb_get_particle_by_hash
     * will return the first hit.  2 hashes generated with the reb_tools_hash hash function have a ~1e-9 chance of clashing.
     */
}

