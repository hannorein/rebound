/**
 * Simulationarchive Viewer
 *
 * This example allows you load in a Simulationarchive and visualize it.
 * You can use the keyboard to step through the individual snapshots.
 * It work with the web-based visualization as well as with OpenGL.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

struct reb_simulationarchive* sa;
void heartbeat(struct reb_simulation* const r);

int64_t current_snapshot = 0;

int key_callback(struct reb_simulation* r, int key){
    switch (key){
        case 262: // right arrow
            current_snapshot++;
            break;
        case 263: // left arrow
            current_snapshot--;
            break;
        case 268: // home
            current_snapshot = 0;
            break;
        case 269: // end
            current_snapshot = sa->nblobs - 1;
            break;
        case 266: // page up
            current_snapshot -= 10;
            break;
        case 267: // page down
            current_snapshot += 10;
            break;
        default: // unknown key
            return 0; // check default keys
    }

    // Update simulation
    if (current_snapshot < 0){
        current_snapshot = 0;
    }
    if (current_snapshot >= sa->nblobs){
        current_snapshot = sa->nblobs - 1;
    }
    r->status = REB_STATUS_SUCCESS; // will trigger reb_simulation_integrate to exit
    return 1;
}

int main(int argc, char* argv[]) {
    if (argc!=2){
        printf("Usage: rebound simulationarchive.bin\n");
        return 1;
    }
    sa = reb_simulationarchive_create_from_file(argv[1]);
    if (!sa){
        printf("Error loading Simulationarchive from file `%s`.\n",argv[1]);
        return 1;
    }

    printf("Simulationarchive loaded from file `%s`.\n",argv[1]);
    printf("Number of snapshots: %lld.\n", sa->nblobs);
    printf("You can step through the Simulationarchive using the following keys in the visualization window:\n");
    printf(" Right arrow: jump to next snapshot\n");
    printf(" Left arrow:  jump to previous snapshot\n");
    printf(" Page down:   jump 10 snapshots foward\n");
    printf(" Page up:     jump 10 snapshots backward\n");
    printf(" Home key:    jump to first snapshot\n");
    printf(" End key:     jump to last snapshot\n\n");

    while(1){
        printf("Loading snapshot %lld.\n", current_snapshot);
        struct reb_simulation* r = reb_simulation_create_from_simulationarchive(sa, current_snapshot);
        if (!r){
            printf("Error loading Simulation from Simulationarchive.\n");
            return 1;
        }

        r->key_callback = key_callback;
        r->status = REB_STATUS_PAUSED;

        // This allows you to connect to the simulation using
        // a web browser by pointing it to http://localhost:1234
        reb_simulation_start_server(r, 1234);

        // Not actually integrating because simulation is paused. 
        reb_simulation_integrate(r, INFINITY);

        if (r->status > 0){ // quit
            reb_simulation_free(r);
            break;
        }

        reb_simulation_free(r);
    }
    reb_simulationarchive_free(sa);
    return 0; 
}
