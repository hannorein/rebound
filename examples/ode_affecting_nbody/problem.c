/**
 * ODE affecting N-body simulation
 * 
 * This examples shows how to integrate arbitrary ODEs
 * and add a back-reaction to N-body particles via an
 * additional_forces function.
 *
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

const double k = 1.; // Constants for the Harmonic Oscillator
const double m = 1.;

struct reb_ode* ho;

void ode_derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    const double omega = sqrt(k/m);
    struct reb_orbit o = reb_orbit_from_particle(ode->r->G, ode->r->particles[1], ode->r->particles[0]);
    double forcing = sin(o.f);
    yDot[0] = y[1]; 
    yDot[1] = -omega*omega*y[0] + forcing;
}

void additional_forces(struct reb_simulation* r){
    struct reb_particle* const particles = r->particles;
    const double coupling = 1e-4;
    // The harmonic oscillator is forcing the planet in the x direction.
    // Note: We are using the variable y1 in the ODE sruct. 
    //       This is the current state of the ODE during the timestep.
    //       The variable y is only updated at the end of a successful timestep.
    particles[1].ax += coupling*ho->y1[0]; 
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();

    reb_simulation_add_fmt(r, "m", 1.);                // Central object
    reb_simulation_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_simulation_move_to_com(r);

    r->integrator = REB_INTEGRATOR_BS;  // Bulirsch-Stoer integrator
    r->ri_bs.eps_rel = 1e-8;            // Relative tolerance
    r->ri_bs.eps_abs = 1e-8;            // Absolute tolerance
    r->dt = 1e-2;

    ho = reb_ode_create(r,2);   // Add an ODE with 2 dimensions
    ho->derivatives = ode_derivatives;          // Right hand side of the ODE
    ho->y[0] = 1;                               // Initial conditions
    ho->y[1] = 0;

    r->additional_forces = additional_forces;

    while(r->t<10){
        reb_simulation_integrate(r, r->t + 0.3);
        printf("y(%.5f) \t = %.5f \n",r->t, ho->y[0]);
    }
    

    reb_ode_free(ho);
    reb_simulation_free(r);

}

