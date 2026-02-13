/**
 * Custom integrator with data storage
 * 
 * This example shows how to use a custom, user-provided integrator
 * that requires internal storage with REBOUND. Here, we implement 
 * the 2nd order leap-frog integrator which combines the first and
 * last Drift steps.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct leapfrog_data {
    int is_synchronized;
};

void leapfrog_step(struct reb_simulation* r){
    if (r->ri_custom.data_size==0){
        // Need to allocate memory on first call
        size_t size = sizeof(struct leapfrog_data);
        r->ri_custom.data = malloc(size);
        r->ri_custom.data_size = size; // size of data in bytes
        // Populate default value
        struct leapfrog_data* data = r->ri_custom.data;
        data->is_synchronized = 1; 
    }

    struct reb_particle* restrict const particles = r->particles;
    struct leapfrog_data* data = r->ri_custom.data;

    // Drift
    if (data->is_synchronized){
        // Half drift step
        for (unsigned int i=0;i<r->N;i++){
            particles[i].x  += (0.5*r->dt) * particles[i].vx;
            particles[i].y  += (0.5*r->dt) * particles[i].vy;
            particles[i].z  += (0.5*r->dt) * particles[i].vz;
        }
        r->t += 0.5*r->dt; // Advance time
    }else{
        // Full drift step
        for (unsigned int i=0;i<r->N;i++){
            particles[i].x  += r->dt * particles[i].vx;
            particles[i].y  += r->dt * particles[i].vy;
            particles[i].z  += r->dt * particles[i].vz;
        }
        r->t += r->dt; // Advance time
    }
    // Kick
    reb_simulation_update_acceleration(r);
    for (unsigned int i=0;i<r->N;i++){
        particles[i].vx += r->dt * particles[i].ax;
        particles[i].vy += r->dt * particles[i].ay;
        particles[i].vz += r->dt * particles[i].az;
    }
    // Skipping second Drift step. Particles are unsynchronized.
    data->is_synchronized = 0; 
}

void leapfrog_synchronize(struct reb_simulation* r){
    if (r->ri_custom.data_size!=0){
        struct leapfrog_data* data = r->ri_custom.data;
        if (data->is_synchronized==0){
            // Drift half step
            struct reb_particle* restrict const particles = r->particles;
            for (unsigned int i=0;i<r->N;i++){
                particles[i].x  += (0.5*r->dt) * particles[i].vx;
                particles[i].y  += (0.5*r->dt) * particles[i].vy;
                particles[i].z  += (0.5*r->dt) * particles[i].vz;
            }
            r->t += 0.5*r->dt; // Advance time
            data->is_synchronized = 1;
        }
    }
}

void leapfrog_reset(struct reb_simulation* r){
    if (r->ri_custom.data_size!=0){
        free(r->ri_custom.data);
        r->ri_custom.data = NULL;
        r->ri_custom.data_size = 0;
    }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_add_fmt(r, "m", 1.);                // Central object
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_simulation_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle 

    // Choose CUSTOM integrator and setup function pointers.
    r->integrator = REB_INTEGRATOR_CUSTOM;
    r->ri_custom.step = leapfrog_step;
    r->ri_custom.synchronize = leapfrog_synchronize;
    r->ri_custom.reset = leapfrog_reset;
    
    // Integrate 5 time units
    // This will automatically call synchronize() at the end.
    r->exact_finish_time = 0; // Allow overshoot. We do not want to change the timestep in last step.
    reb_simulation_integrate(r, 5.0);
    printf("is_synchronized = %d\n", ((struct leapfrog_data*)r->ri_custom.data)->is_synchronized);
    
    // Step forward 10 steps
    reb_simulation_steps(r, 10);
    // Simulation is not synchronized afterwards.
    printf("is_synchronized = %d\n", ((struct leapfrog_data*)r->ri_custom.data)->is_synchronized);
    // Synchronize manually
    reb_simulation_synchronize(r);
    printf("is_synchronized = %d\n", ((struct leapfrog_data*)r->ri_custom.data)->is_synchronized);

    // Step forward 10 more steps
    reb_simulation_steps(r, 10);
    // We can save the simulation to a file in the unsynchronized state.
    // This automatically saves everything stored in r->ri_custom.data to the file.
    reb_simulation_save_to_file(r, "out.bin");
    // Free the original simulation. This will call leapfrog_reset() which is responsible for freeing r->ri_custom.data.
    reb_simulation_free(r);

    // Load the simulation back from the file.
    r = reb_simulation_create_from_file("out.bin", -1);
    // Need to set the function pointers manually. 
    // Data is automatically restored.
    r->ri_custom.step = leapfrog_step;
    r->ri_custom.synchronize = leapfrog_synchronize;
    r->ri_custom.reset = leapfrog_reset;
    // Synchronize manually
    reb_simulation_synchronize(r);
    printf("is_synchronized = %d\n", ((struct leapfrog_data*)r->ri_custom.data)->is_synchronized);

    // cleanup
    reb_simulation_free(r);
}

