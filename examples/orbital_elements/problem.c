/**
 * Orbital Elements
 *
 * This is a simple overview of the helper functions in REBOUND
 * for using orbital elements.  For an in-depth discussion and
 * examples, see ipython_examples/OrbitalElements.ipynb.  The
 * C function calls are more explicit (see below), but the numerical
 * issues and conventions are the same in Python and C.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    struct reb_particle p; 
    p.m      = 1.;    
    reb_simulation_add(r, p); 

    // Adding a particle with orbital elements requires the following 7 things
    // There's more flexibility in python for passing different orbital elements
    // for edge cases.  The user has to calculate these manually in C and pass
    // the elements below
    
    struct reb_particle primary = r->particles[0];
    double m = 0.;
    double a = 0.1;
    double e = 0.2;
    double inc = 0.3;
    double Omega = 0.4;
    double omega = 0.5;
    double f = 0.6;

    struct reb_particle p2 = reb_particle_from_orbit(r->G, primary, m, a, e, inc, Omega, omega, f);
    reb_simulation_add(r,p2);

    struct reb_orbit o= reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);
    
    printf("a = %.16e\n", o.a);
    printf("e = %.16e\n", o.e);
    printf("inc = %.16e\n", o.inc);
    printf("Omega = %.16e\n", o.Omega);
    printf("omega = %.16e\n", o.omega);
    printf("f = %.16e\n", o.f);

    // There are also versions of the two functions above that let you pass an error integer pointer
    // to diagnose / catch errors.  You can find error codes in the documentation for the functions
    // at https://rebound.readthedocs.org/
    
    int err = 0;
    e = 1.001;

    p2 = reb_particle_from_orbit_err(r->G, primary, m, a, e, inc, Omega, omega, f, &err);
    if(err == 3){            // error code for bound orbit with e > 1
        e = 1.-1.e-15;        // set to just less than 1
        p2 = reb_particle_from_orbit_err(r->G, primary, m, a, e, inc, Omega, omega, f, &err);
    }
    reb_simulation_add(r,p2);
    
    o= reb_orbit_from_particle(r->G, r->particles[2], r->particles[0]);

    printf("\n\na = %.16e\n", o.a);
    printf("e = %.16e\n", o.e);
    printf("inc = %.16e\n", o.inc);
    printf("Omega = %.16e\n", o.Omega);
    printf("omega = %.16e\n", o.omega);
    printf("f = %.16e\n", o.f);

    reb_simulation_free(r);
}
