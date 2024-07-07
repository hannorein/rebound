/**
 * Screenshots
 *
 * This example shows how to take a screenshot of a REBOUND 
 * simulation at the beginning of the simulation and once 
 * every 400 timesteps during the integration. You need to
 * compile REBOUND with SERVER=1 and connect a webbrowser
 * to the simulation. Only then can you take screenshots. 
 * You might also be interested in the examples:
 * 1) animation_saturns_rings and
 * 2) animation_solar_system.
 * They show how one can programatically change the
 * visualization. You can combine this with taking
 * screenshots if you want to record animations of REBOUND
 * simulations.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* const r){
    if (r->steps_done%400==0){ // Every 400 timesteps
        int id = r->steps_done/400;
        char filename[1024];
        sprintf(filename, "screenshot_%05d.png", id);
        if (reb_simulation_output_screenshot(r, filename)){
            printf("Screenshot saved: %s\n", filename);
        }
    }
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->heartbeat = heartbeat;

    // Initial conditions
    reb_simulation_add_fmt(r, "m", 1.); // star
    reb_simulation_add_fmt(r, "m a", 1e-3, 1.); // planet 1
    reb_simulation_add_fmt(r, "m a", 1e-3, 2.); // planet 2
    reb_simulation_move_to_com(r);
   
    // Start the web server. Make sure you point your 
    // webbrowser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
    
    // Manually take a screenshot of the initial conditions.
    // The program will pause here until you connect a
    // webbrowser and a screenshot can be taken.
    reb_simulation_output_screenshot(r, "screenshot_initial.png");
    printf("Screenshot saved: screenshot_initial.png\n");

    // Start integration.
    reb_simulation_integrate(r, 10);
    
    // Manually take a screenshot of the final simulation
    reb_simulation_output_screenshot(r, "screenshot_final.png");
    printf("Screenshot saved: screenshot_final.png\n");

    // Cleanup
    reb_simulation_free(r);
}

