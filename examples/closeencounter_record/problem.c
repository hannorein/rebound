/**
 * Detect and record close encounters
 *
 * This example integrates a densely packed planetary system 
 * which becomes unstable on a timescale of only a few orbits. 
 * The example is identical to the `close_encounter` sample, except that 
 * the collisions are recorded and written to a file. What kind of collisions
 * are recorded can be easily modified. It is also possible to implement some
 * additional physics whenever a collision has been detection (e.g. fragmentation).
 * The collision search is by default a direct search, i.e. O(N^2) but can be
 * changed to a tree by using the `collisions_tree.c` module.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

// Define our own collision resolve function, which will only record collisions but not change any of the particles.        
int collision_record_only(struct reb_simulation* const r, struct reb_collision c){
    double delta_t = 2.*M_PI;     
    struct reb_particle* particles = r->particles;
    const double t = r->t;

    // only record a maximum of one collision per year per particle
    if ( particles[c.p1].lastcollision+delta_t < t  &&  particles[c.p2].lastcollision+delta_t < t ){
        particles[c.p1].lastcollision = t; 
        particles[c.p2].lastcollision = t;
        printf("\nCollision detected.\n");        
        FILE* of = fopen("collisions.txt","a+");        // open file for collision output
        fprintf(of, "%e\t", t);                    // time
        fprintf(of, "%e\t", (particles[c.p1].x+particles[c.p2].x)/2.);    // x position
        fprintf(of, "%e\t", (particles[c.p1].y+particles[c.p2].y)/2.);    // y position
        fprintf(of, "\n");
        fclose(of);                        // close file
    }
    return 0;
}


void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10.*2.*M_PI)){  
        reb_output_timing(r, 0);
    }
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    r->dt             = 0.1*2.*M_PI;                // initial timestep
    r->integrator         = REB_INTEGRATOR_IAS15;
    r->collision        = REB_COLLISION_DIRECT;
    r->collision_resolve     = collision_record_only;        // Set function pointer for collision recording.
    r->heartbeat        = heartbeat;
    r->usleep        = 10000;                // Slow down integration (for visualization only)

    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0;                            // Star is pointmass
    reb_add(r, star);
    
    // Add planets
    int N_planets = 7;
    for (int i=0; i<N_planets; i++){
        double a = 1.+(double)i/(double)(N_planets-1);        // semi major axis
        double v = sqrt(1./a);                     // velocity (circular orbit)
        struct reb_particle planet = {0};
        planet.m = 1e-4; 
        double rhill = a * pow(planet.m/(3.*star.m),1./3.);    // Hill radius
        planet.r = rhill;                    // Set planet radius to hill radius 
                                    // A collision is recorded when planets get within their hill radius
                                    // The hill radius of the particles might change, so it should be recalculated after a while
        planet.lastcollision = 0; 
        planet.x = a;
        planet.vy = v;
        reb_add(r, planet); 
    }
    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.

    reb_integrate(r, INFINITY);
}
