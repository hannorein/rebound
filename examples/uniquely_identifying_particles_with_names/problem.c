/**
 * How to use names to identify particles
 *
 * This example shows how to assign names to particles
 * and how to access particles using names.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char* argv[]){
    /* Several events can make particles move around in memory. This means the user should not assume they can always refer to 
     * the same particle through r->particles[index] or through a pointer to the particle that is set at the beginning of the 
     * simulation.  The reliable way to access particles is through names.*/

    struct reb_simulation* r = reb_simulation_create();

    /* Let's add our first particle with the name "Sun" */
    struct reb_particle p = {0};
    p.m = 1.;
    p.name = "Sun";
    reb_simulation_add(r, p);   // The string "Sun" gets copied by the simulation when the particle is added.

    struct reb_particle* sun = reb_simulation_get_particle_by_name(r, "Sun");
    if (sun){ // Particle found?
        printf("We can now reference the particle like this: m=%f, or like this: m=%f\n", r->particles[0].m, sun->m);
        printf("We can also print out the particle's name: %s\n", sun->name);
    }

    // We now add 200 particles. Each has its own name.
    for (int i=1; i <= 200; i++){
        struct reb_particle tp = {0};
        tp.x = i;       // put particles progressively farther away
        char name[256];
        sprintf(name, "testparticle %d", i);
        tp.name = name;
        reb_simulation_add(r, tp);  // Name gets copied by the simulation. Memory is managed by the simulation.
    }

    /* The advantage of names is that if particles are ejected or otherwise removed from the simulation (or you are using 
     * the tree code), the indices in the particles array will get scrambled.  Accessing particles through their name 
     * guarantees you get back the particle you intended.
     */

    /* Initially, the particle with index 200 is "testparticle 200" */

    printf("Initially, the particle with index 200 is: r->particles[200].name\"=%s\"\n", r->particles[200].name);

    reb_simulation_remove_particle_by_name(r, "Sun");

    /* The remove function has moved particles[200] to index 0:*/

    printf("r->particles[0].name=\"%s\"\n", r->particles[0].name);

    /* Rather than worry about the internals of what the remove function, the tree code etc. do, we can get it by name. 
     * When you are not sure, assigning particles names and accessing them through them is always safe. */

    struct reb_particle* p200 = reb_simulation_get_particle_by_name(r, "testparticle 200");

    printf("p200->name=\"%s\"\n", p200->name); 

    /* Note that if the particle is not found in the simulation, reb_simulation_get_particle_by_name returns a NULL pointer.
     * This allows you to check if particles are still in the simulation when you don't know ahead of time, but means that if you 
     * mistakenly access a removed particle, you'll get a segmentation fault.*/

    sun = reb_simulation_get_particle_by_name(r, "Sun");

    if (sun == NULL){
        printf("Whoops!  Already removed particle.\n");
    }
    else{
        printf("Mass of \"Sun\": %f\n\n", sun->m); // would cause segmentation fault in this case
    }

    /* If two particles share the same name, reb_simulation_get_particle_by_name returns a pointer to either one. 
     * Which one is undefined.
     */

    reb_simulation_free(r);
}

