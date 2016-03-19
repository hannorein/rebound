/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
    struct reb_simulation* sim = reb_create_simulation();

    struct reb_particle p = {0.};
    p.m = 1;
    reb_add(sim,p);
    p = reb_tools_orbit_to_particle(1, sim->particles[0],0.,
            1., //a
            0.2, //e
            0.4, // i
            0.4, // Omega
            1.4, // omega
            3.); //f
    reb_add(sim,p);

    double a, lambda, k, h, ix, iy;

    reb_integrate(sim,10.);
   
    p = sim->particles[1]; 
    reb_particle_to_pal(1.,sim->particles[0],p,&a,&lambda,&k,&h,&ix,&iy);
    struct reb_particle pp = reb_pal_to_particle(1.,sim->particles[0],0.,a,lambda,k,h,ix,iy);
    printf("%f %f %f %f %f %f\n",p.x,p.y,p.z, p.vx,p.vy,p.vz);
    printf("%f %f %f %f %f %f\n",pp.x,pp.y,pp.z, pp.vx,pp.vy,pp.vz);

}

