/**
 * Example usage of reb_simulation_add_fmt()
 * 
 * The reb_simulation_add_fmt() function can be used to add a particle 
 * to the simulation by specifying the particle's coordinates
 * in a variety of formats. 
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    // Central object at origin with mass=1
    reb_simulation_add_fmt(r, "m", 1.); 
    
    // Massless particle on circular orbit with a=1
    reb_simulation_add_fmt(r, "a", 1.); 

    // A Jupiter mass planet with a=2 and e=0.1
    reb_simulation_add_fmt(r, "m a e", 1e-3, 2., 0.1); 

    // By default Jacobi coordinates are used, here we use
    // heliocentric coordinates by specifying the primary.
    reb_simulation_add_fmt(r, "a e primary", 0.1, 0.3, r->particles[0]); 

    // The function supports any arbitrary combination of 
    // orbital parameters as long as it's physically meaningful. 
    reb_simulation_add_fmt(r, "a e omega", 3., 0.1, M_PI);
    reb_simulation_add_fmt(r, "a e pomega", 4., 0.1, M_PI/2.);
    reb_simulation_add_fmt(r, "P h k", 365.25, 0.01, 0.02);
    
    // Non-physical parameter combinations raise errors 
    reb_simulation_add_fmt(r, "a e", -1., 0.1);
    reb_simulation_add_fmt(r, "a e", -1., 0.1);

    // Cartesian coordinates are supported as well
    reb_simulation_add_fmt(r, "m x vy", 1e-3, 9., 0.3);

    /** Supported parameters are: 
     *       m:          Mass  (Default: 0)
     *       x, y, z:    Positions in Cartesian coordinates  (Default: 0)
     *       vx, vy, vz: Velocities in Cartesian coordinates (Default: 0)
     *       primary:    Primary body for converting orbital elements to cartesian (Default: center of mass of the particles in the passed simulation, i.e., this will yield Jacobi coordinates as one progressively adds particles) 
     *       a:          Semimajor axis (a or P required if passing orbital elements)
     *       P:          Orbital period (a or P required if passing orbital elements)
     *       e:          Eccentricity                (Default: 0)
     *       inc:        Inclination                 (Default: 0)
     *       Omega:      Longitude of ascending node (Default: 0)
     *       omega:      Argument of pericenter      (Default: 0)
     *       pomega:     Longitude of pericenter     (Default: 0)
     *       f:          True anomaly                (Default: 0)
     *       M:          Mean anomaly                (Default: 0)
     *       E:          Eccentric anomaly           (Default: 0)
     *       l:          Mean longitude              (Default: 0)
     *       theta:      True longitude              (Default: 0)
     *       T:          Time of pericenter passage  
     *       h:          h variable, see Pal (2009) for a definition  (Default: 0)
     *       k:          k variable, see Pal (2009) for a definition  (Default: 0)
     *       ix:         ix variable, see Pal (2009) for a definition  (Default: 0)
     *       iy:         iy variable, see Pal (2009) for a definition  (Default: 0)
     *       r:          Particle radius
     *
     * Note that it is important to pass floating point numbers and not integers to reb_simulation_add_fmt(). 
     * For example: 
     *     reb_simulation_add_fmt(r, "a", 1) 
     * will lead to undefined behaviour. Instead, use 
     *     reb_simulation_add_fmt(r, "a", 1.0) 
     */

    // Run a test simulation
    reb_simulation_move_to_com(r);
    reb_simulation_integrate(r,100.);
}

