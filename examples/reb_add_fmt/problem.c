/**
 * Example usage of reb_add_fmt()
 * 
 * The reb_add_fmt() function can be used to add particles 
 * to the simulation by specifying the particles coordinates
 * in a variety of formats. 
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    
    // Central object at origin with mass=1
    reb_add_fmt(r, "m", 1.); 
    
    // Massless particle on circular orbit with a=1
    reb_add_fmt(r, "a", 1.); 

    // A Jupiter mass planet with a=2 and e=0.1
    reb_add_fmt(r, "m a e", 1e-3, 2., 0.1); 

    // By default Jacobi coordinates are used, here we use
    // heliocentric coordinates by specifying the primary.
    reb_add_fmt(r, "a e primary", 0.1, 0.3, r->particles[0]); // By default 

    // The function supports any arbitrary combination of 
    // orbital parameters as long as it's physically meaningful. 
    reb_add_fmt(r, "a e omega", 3., 0.1, M_PI);
    reb_add_fmt(r, "a e pomega", 4., 0.1, M_PI/2.);
    reb_add_fmt(r, "P h k", 365.25, 0.01, 0.02);
    
    // Non-physical parameter combinations raise errors 
    reb_add_fmt(r, "a e", -1., 0.1);
    reb_add_fmt(r, "a e", -1., 0.1);

    // Cartesian coordinates are supported as well
    reb_add_fmt(r, "m x vy", 1e-3, 9., 0.3);
 
    /** Note that it is important to pass floating point numbers, not 
     * integers, to reb_add_fmt(). For example: 
     *     reb_add_fmt(r, "a", 1) 
     * will lead to undefined behaviour. Use 
     *     reb_add_fmt(r, "a", 1.0) 
     * instead.
     */

    // Run a test simulation
    reb_move_to_com(r);
    reb_integrate(r,100.);
}

