/**
 * Custom integrator with state storage
 * 
 * This example shows how to use a custom, user-provided integrator
 * that requires internal storage with REBOUND. Here, we implement 
 * the 2nd order leap-frog integrator which combines the first and
 * last Drift steps. We also record some diagnostic data in the 
 * integrator state to demonstrate how to build a custom integrator
 * which requires its own datastorage.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct custom_integrator_state {
    // Here we just log the number of steps
    // but this could contain any data or array 
    // as required by the integrator.
    int number_of_steps;
    int number_of_synchronizations;
};

// By defining the layout of the state in the form of a reb_binary_data_field_descriptor,
// REBOUND can automatically safe and restore integrator state to a Simulationarchive file.
const struct reb_binarydata_field_descriptor leapfrog_field_descriptor_list[] = {
    { 11042, REB_INT, "number_of_steps",            offsetof(struct custom_integrator_state, number_of_steps), 0, 0, 0},
    { 11043, REB_INT, "number_of_synchronizations", offsetof(struct custom_integrator_state, number_of_synchronizations), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* leapfrog_create(){
    struct custom_integrator_state* leapfrog = malloc(sizeof(struct custom_integrator_state));
    // Default values 
    leapfrog->number_of_steps = 0;  
    leapfrog->number_of_synchronizations = 0;
    return leapfrog;
}

void leapfrog_free(void* p){
    struct custom_integrator_state* leapfrog = p;
    // Free any data allocated by integrator
    // Here, it is just the state itself.
    free(leapfrog);
}

void leapfrog_step(struct reb_simulation* r, void* p){
    struct reb_particle* restrict const particles = r->particles;
    struct custom_integrator_state* leapfrog = p;

    // Drift
    if (r->is_synchronized){
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
    r->is_synchronized = 0; 

    leapfrog->number_of_steps++;
}

void leapfrog_synchronize(struct reb_simulation* r, void* p){
    if (r->is_synchronized==0){
        struct custom_integrator_state* leapfrog = p;
        // Drift half step
        struct reb_particle* restrict const particles = r->particles;
        for (unsigned int i=0;i<r->N;i++){
            particles[i].x  += (0.5*r->dt) * particles[i].vx;
            particles[i].y  += (0.5*r->dt) * particles[i].vy;
            particles[i].z  += (0.5*r->dt) * particles[i].vz;
        }
        r->t += 0.5*r->dt; // Advance time
        r->is_synchronized = 1;
        leapfrog->number_of_synchronizations++;
    }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    reb_simulation_add_fmt(r, "m", 1.);                // Central object
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_simulation_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle 

    // Choose CUSTOM integrator and setup function pointers.
    struct reb_integrator custom = {
        .step = leapfrog_step,
        .synchronize = leapfrog_synchronize,
        .create = leapfrog_create,
        .free = leapfrog_free,
        .field_descriptor_list = leapfrog_field_descriptor_list,
    };
    reb_integrator_register(custom, "custom-integrator");
    reb_simulation_set_integrator(r, "custom-integrator");
    
    // Integrate 5 time units
    // This will automatically call synchronize() at the end.
    r->exact_finish_time = 0; // Allow overshoot. We do not want to change the timestep in last step.
    reb_simulation_integrate(r, 5.0);
    struct custom_integrator_state* leapfrog = r->integrator.state;
    printf("number_of_steps = %d\n", leapfrog->number_of_steps);
    printf("number_of_synchronizations = %d\n", leapfrog->number_of_synchronizations);
    
    // Step forward 10 steps
    reb_simulation_steps(r, 10);
    // Simulation is not synchronized right now.
    printf("is_synchronized = %d\n", r->is_synchronized);
    // Synchronize manually
    reb_simulation_synchronize(r);
    printf("is_synchronized = %d\n", r->is_synchronized);

    // Step forward 10 more steps
    reb_simulation_steps(r, 10);

    // We next save the simulation to a file.
    // This automatically saves everything stored in r->integrator.state to the file.
    // Note that it is in the unsynchronized state.
    reb_simulation_save_to_file(r, "out.bin");
    // Free the original simulation. This will call leapfrog_reset() which is responsible for freeing r->integrator.state.
    reb_simulation_free(r);

    // Load the simulation back from the file.
    // Custom integrator needs to be registered before loading file. 
    // Integrator state is restored automatically if field_descriptor_list is provided.
    r = reb_simulation_create_from_file("out.bin", -1);
    leapfrog = r->integrator.state;
    printf("number_of_steps = %d\n", leapfrog->number_of_steps);
    printf("number_of_synchronizations = %d\n", leapfrog->number_of_synchronizations);

    
    // Simulation is not synchronized right now.
    printf("is_synchronized = %d\n", r->is_synchronized);
    // Synchronize manually
    reb_simulation_synchronize(r);
    printf("is_synchronized = %d\n", r->is_synchronized);

    // cleanup
    reb_simulation_free(r);
}

