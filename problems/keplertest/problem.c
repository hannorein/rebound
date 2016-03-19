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
    p.m = 0.;
    p.x = 1.;
    p.vy = 1.;
    reb_add(sim,p);

    double a, lambda, k, h, ix, iy;
    reb_particle_to_pal(1.,sim->particles[0],sim->particles[1],&a,&lambda,&k,&h,&ix,&iy);
    printf("%f %f %f %f %f %f\n",a, lambda, k, h, ix, iy);

}

