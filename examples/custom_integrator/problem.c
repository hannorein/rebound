/**
 * Custom integrator
 * 
 * This example shows how to use a custom, user-provided integrator
 * with REBOUND. Here, we implement the 2nd order leap-frog integrator.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

void leapfrog_step(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    // Drift
    for (unsigned int i=0;i<r->N;i++){
        particles[i].x  += (0.5*r->dt) * particles[i].vx;
        particles[i].y  += (0.5*r->dt) * particles[i].vy;
        particles[i].z  += (0.5*r->dt) * particles[i].vz;
    }
    r->t += 0.5*r->dt; // Advance time
    // Kick
    reb_simulation_update_acceleration(r);
    for (unsigned int i=0;i<r->N;i++){
        particles[i].vx += r->dt * particles[i].ax;
        particles[i].vy += r->dt * particles[i].ay;
        particles[i].vz += r->dt * particles[i].az;
    }
    // Drift
    for (unsigned int i=0;i<r->N;i++){
        particles[i].x  += (0.5*r->dt) * particles[i].vx;
        particles[i].y  += (0.5*r->dt) * particles[i].vy;
        particles[i].z  += (0.5*r->dt) * particles[i].vz;
    }
    r->t += 0.5*r->dt; // Advance time
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_add_fmt(r, "m", 1.);                // Central object
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_simulation_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle 

    // Create a copy of the simulation and use the built-in leapfrog integrator
    // as a comparison.
    struct reb_simulation* r_copy= reb_simulation_copy(r);
    r_copy->integrator = REB_INTEGRATOR_LEAPFROG;

    // Choose CUSTOM integrator and setup function pointers.
    r->integrator = REB_INTEGRATOR_CUSTOM;
    r->ri_custom.step = leapfrog_step;
    
    // Integrate both simulations
    reb_simulation_integrate(r,10.);
    reb_simulation_integrate(r_copy,10.);

    // Compare the final coordinates of both simulations.
    for (int i=0; i<r->N; i++){
        struct reb_particle p1 = r->particles[i];
        struct reb_particle p2 = r_copy->particles[i];
        // They are be identical.
        assert(p1.x==p2.x);
        assert(p1.y==p2.y);
        assert(p1.z==p2.z);
    }

    printf("Simulations complete. Tests passed.\n");

    // Cleanup 
    reb_simulation_free(r);
    reb_simulation_free(r_copy);
}

