/**
 * Variational Equations
 * 
 * This example shows how to use first and second
 * order variational equations.
 * See also https://github.com/hannorein/rebound/blob/master/ipython_examples/VariationalEquations.ipynb and Rein and Tamayo (2016).
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

// This function creates a simulation with one star, one planet and one test particle.
struct reb_simulation* create_sim(){
    struct reb_simulation* r = reb_simulation_create();
    r->integrator = REB_INTEGRATOR_IAS15; // First and second order  variational equations supported in IAS15.
    // r->integrator = REB_INTEGRATOR_BS; // First and second order  variational equations supported in BS.
    // r->integrator = REB_INTEGRATOR_WHFAST;  Only first order variational equations supported in WHFast.
    struct reb_particle star = {0.};
    star.m = 1;
    reb_simulation_add(r, star);
    struct reb_particle planet = reb_particle_from_orbit(1.,star,1e-3,1.,0.,0.,0.,0.,0.);
    reb_simulation_add(r, planet);
    struct reb_particle testparticle = reb_particle_from_orbit(1.,star,0.,1.7,0.1,0.2,0.3,0.4,0.5);
    reb_simulation_add(r, testparticle);
    return r;
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r;
    int var_i, var_ii;

    // We first integrate the vanilla simulation forward in time and look at the position of the testparticle at the end of the simulation.
    r = create_sim();
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100:                             %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_simulation_free(r);
    
    // Next, we shift the planet's initial x coordinate by a small amount and integrate the system again up til t=100. 
    double DeltaX = 0.001;
    printf("\nShifting planet's x coordinate by %f.\n", DeltaX);
    r = create_sim();
    r->particles[1].x += DeltaX;
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 in shifted simulation:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_simulation_free(r);
    
    // Instead of shifting the initial x coordinate, we can also use variational equations. 
    r = create_sim();
    var_i = reb_simulation_add_variation_1st_order(r, -1);   // The -1 means we vary a particle which is not a testparticle  and therefore can influence other particles
    
    // By default all components of variational particles are initialized to zero.
    // We are interested in shifting the planet's x coordinates and thus initialize the x coordinate of the variational particle to 1.
    r->particles[var_i+1].x = 1.;           
    reb_simulation_integrate(r,100.);
    // After the integration ran, we can estimate where the test particle would have been had we shifted the inner planet's initial x coordinate.
    printf("Position of testparticle at t=100 using 1st order var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i+2].x,r->particles[2].y+DeltaX*r->particles[var_i+2].y);
    reb_simulation_free(r);

    // Better yet, we can use second order variational particles.
    r = create_sim();
    var_i = reb_simulation_add_variation_1st_order(r, -1);
    var_ii = reb_simulation_add_variation_2nd_order(r, -1, var_i, var_i);
    r->particles[var_i+1].x = 1.;
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd order var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i+2].x+DeltaX*DeltaX/2.*r->particles[var_ii+2].x,r->particles[2].y+DeltaX*r->particles[var_i+2].y+DeltaX*DeltaX/2.*r->particles[var_ii+2].y);
    reb_simulation_free(r);

    
    // We now do the same as above, but vary the testparticle's position 
    printf("\nShifting testparticle's x coordinate by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2].x += DeltaX;
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 in shifted simulation:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_simulation_free(r);
    
    r = create_sim();
    var_i = reb_simulation_add_variation_1st_order(r, 2); // The 2 corresponds to the index of the testparticle that we vary.
    r->particles[var_i].x = 1.;
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st order var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_simulation_free(r);

    r = create_sim();
    var_i = reb_simulation_add_variation_1st_order(r, 2);
    var_ii = reb_simulation_add_variation_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i].x = 1.;
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd order var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_simulation_free(r);

    
    // Instead of varying cartesian coordinates, we can also vary orbital elements. 
    printf("\nShifting planet's semi-major axis by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2] = reb_particle_from_orbit(1.,r->particles[0],0.,1.7+DeltaX,0.1,0.2,0.3,0.4,0.5);
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 in shifted simulation:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_simulation_free(r);

    r = create_sim();
    var_i = reb_simulation_add_variation_1st_order(r, 2);
    // The function that sets up the variational particle gets the same orbital parameters as the original particle.
    r->particles[var_i] = reb_particle_derivative_a(1.,r->particles[0],r->particles[2]);
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st order var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_simulation_free(r);
    
    r = create_sim();
    var_i = reb_simulation_add_variation_1st_order(r, 2);
    var_ii = reb_simulation_add_variation_2nd_order(r, 2, var_i, var_i);
    // first derivative with respect to a
    r->particles[var_i] = reb_particle_derivative_a(1.,r->particles[0],r->particles[2]);
    // second derivative with respect to a
    r->particles[var_ii] = reb_particle_derivative_a_a(1.,r->particles[0],r->particles[2]);
    reb_simulation_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd order var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_simulation_free(r);


}
