/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include "output.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m", 1.);                // Central object
    reb_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_add_fmt(r, "a e", 1.4, 0.1);        // Massless test particle 

    reb_steps(r,1);

    char* buf = NULL;
    size_t size = 0;
    reb_output_binary_to_stream(r, &buf, &size);
    


    struct reb_simulationarchive* sa = malloc(sizeof(struct reb_simulationarchive));
    reb_read_simulationarchive_from_buffer_with_messages(sa, buf, size, NULL, NULL);

    
    struct reb_simulation* r2 = reb_create_simulation();
    enum reb_input_binary_messages w =0;
    reb_create_simulation_from_simulationarchive_with_messages(r2,sa,-1,&w);

    reb_close_simulationarchive(sa);
    free(buf);
    printf("sizeof(size_t)= %ld\n", sizeof(unsigned int));

    reb_free_simulation(r);
    reb_free_simulation(r2);
}

