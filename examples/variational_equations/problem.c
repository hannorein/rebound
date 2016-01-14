/**
 * Variational Equations
 * 
 * This example shows how to use first and second
 * order variational equations.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>


struct reb_simulation* create_sim(){
	struct reb_simulation* r = reb_create_simulation();
    // r->integrator = REB_INTEGRATOR_WHFAST;  Only first order variational equations supported in WHFast.
    struct reb_particle star = {0.};
    star.m = 1;
    reb_add(r, star);
    struct reb_particle planet = reb_tools_orbit_to_particle(1.,star,1e-3,1.,0.,0.,0.,0.,0.);
    reb_add(r, planet);
    struct reb_particle testparticle = reb_tools_orbit_to_particle(1.,star,0.,1.7,0.1,0.2,0.3,0.4,0.5);
    reb_add(r, testparticle);
    reb_move_to_com(r);
    return r;
}

int main(int argc, char* argv[]) {
	struct reb_simulation* r;
    double DeltaX = 0.001;
    int var_i, var_ii;

    r = create_sim();
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100:                       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);
    
    
    printf("\nShifting planet's x coordinate by %f.\n", DeltaX);
    r = create_sim();
    r->particles[1].x += DeltaX;
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 in shifted simulation: %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, -1);
    r->particles[var_i+1].x = 1.;
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i+2].x,r->particles[2].y+DeltaX*r->particles[var_i+2].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, -1);
    var_ii = reb_add_var_2nd_order(r, -1, var_i, var_i);
    r->particles[var_i+1].x = 1.;
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i+2].x+DeltaX*DeltaX/2.*r->particles[var_ii+2].x,r->particles[2].y+DeltaX*r->particles[var_i+2].y+DeltaX*DeltaX/2.*r->particles[var_ii+2].y);
    reb_free_simulation(r);

    
    
    printf("\nShifting testparticle's x coordinate by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2].x += DeltaX;
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 in shifted simulation: %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    r->particles[var_i].x = 1.;
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    var_ii = reb_add_var_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i].x = 1.;
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_free_simulation(r);

    
    
    printf("\nShifting planet's a by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2] = reb_tools_orbit_to_particle(1.,r->particles[0],0.,1.7+DeltaX,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using shifted a:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    r->particles[var_i] = reb_tools_orbit_to_particle_da(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    var_ii = reb_add_var_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i] = reb_tools_orbit_to_particle_da(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
    r->particles[var_ii] = reb_tools_orbit_to_particle_dda(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_free_simulation(r);


    
    DeltaX = 0.01; 
    printf("\nShifting planet's e by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2] = reb_tools_orbit_to_particle(1.,r->particles[0],0.,1.7,0.1+DeltaX,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using shifted e:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    r->particles[var_i] = reb_tools_orbit_to_particle_de(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    var_ii = reb_add_var_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i] = reb_tools_orbit_to_particle_de(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
    r->particles[var_ii] = reb_tools_orbit_to_particle_dde(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_free_simulation(r);
   

    DeltaX = 0.05; 
    printf("\nShifting planet's i by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2] = reb_tools_orbit_to_particle(1.,r->particles[0],0.,1.7,0.1,0.2+DeltaX,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using shifted i:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    r->particles[var_i] = reb_tools_orbit_to_particle_di(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    var_ii = reb_add_var_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i] = reb_tools_orbit_to_particle_di(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
    r->particles[var_ii] = reb_tools_orbit_to_particle_ddi(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_free_simulation(r);
   

    DeltaX = 0.05; 
    printf("\nShifting planet's Omega by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2] = reb_tools_orbit_to_particle(1.,r->particles[0],0.,1.7,0.1,0.2,0.3+DeltaX,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using shifted i:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    r->particles[var_i] = reb_tools_orbit_to_particle_dOmega(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    var_ii = reb_add_var_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i] = reb_tools_orbit_to_particle_dOmega(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
    r->particles[var_ii] = reb_tools_orbit_to_particle_ddOmega(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_free_simulation(r);
   

    DeltaX = 0.05; 
    printf("\nShifting planet's omega by %f.\n", DeltaX);
    r = create_sim();
    r->particles[2] = reb_tools_orbit_to_particle(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4+DeltaX,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using shifted i:       %.8f %.8f\n",r->particles[2].x,r->particles[2].y);
    reb_free_simulation(r);

    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    r->particles[var_i] = reb_tools_orbit_to_particle_domega(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 1st var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x,r->particles[2].y+DeltaX*r->particles[var_i].y);
    reb_free_simulation(r);
    
    r = create_sim();
    var_i = reb_add_var_1st_order(r, 2);
    var_ii = reb_add_var_2nd_order(r, 2, var_i, var_i);
    r->particles[var_i] = reb_tools_orbit_to_particle_domega(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
    r->particles[var_ii] = reb_tools_orbit_to_particle_ddomega(1.,r->particles[0],0.,1.7,0.1,0.2,0.3,0.4,0.5);
	reb_integrate(r,100.);
    printf("Position of testparticle at t=100 using 2nd var. eqs.:   %.8f %.8f\n",r->particles[2].x+DeltaX*r->particles[var_i].x+DeltaX*DeltaX/2.*r->particles[var_ii].x,r->particles[2].y+DeltaX*r->particles[var_i].y+DeltaX*DeltaX/2.*r->particles[var_ii].y);
    reb_free_simulation(r);
}

