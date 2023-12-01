/**
 * OpenMP example.
 *
 * A self-gravitating disc is integrated using
 * the leap frog integrator and direct summation.
 * Shared memory parallelization using OpenMP 
 * is enabled in the Makefile.
 *
 * Note that you need a compiler which supports 
 * OpenMP to run this example. By default, the 
 * OSX compilers from Apple do currently not
 * support OpenMP. You can install the GNU 
 * C compilers easily with homebrew. Look at the
 * Makefile of this example to see how you can setup
 * the parameters to compile REBOUND on both OSX
 * and Linux.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "rebound.h"
#include "tools.h"
#include "output.h"


void run_sim(){
    struct reb_simulation* const r = reb_simulation_create();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_LEAPFROG;
    r->gravity    = REB_GRAVITY_BASIC;
    r->boundary    = REB_BOUNDARY_OPEN;
    r->opening_angle2    = 1.5;    // This constant determines the accuracy of the tree code gravity estimate.
    r->G         = 1;        
    r->softening     = 0.02;        // Gravitational softening length
    r->dt         = 3e-2;        // Timestep
    const double boxsize = 10.2;
    reb_simulation_configure_box(r,boxsize,1,1,1);

    // Setup particles
    double disc_mass = 2e-1;    // Total disc mass
    int N = 2000;            // Number of particles
    // Initial conditions
    struct reb_particle star = {0};
    star.m         = 1;
    reb_simulation_add(r, star);
    for (int i=0;i<N;i++){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(r, boxsize/10.,boxsize/2./1.2,-1.5);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt.x         = a*cos(phi);
        pt.y         = a*sin(phi);
        pt.z         = a*reb_random_normal(r, 0.001);
        double mu     = star.m + disc_mass * (pow(a,-3./2.)-pow(boxsize/10.,-3./2.))/(pow(boxsize/2./1.2,-3./2.)-pow(boxsize/10.,-3./2.));
        double vkep     = sqrt(r->G*mu/a);
        pt.vx         =  vkep * sin(phi);
        pt.vy         = -vkep * cos(phi);
        pt.vz         = 0;
        pt.m         = disc_mass/(double)N;
        reb_simulation_add(r, pt);
    }

    reb_simulation_integrate(r, 1.0);
    reb_simulation_free(r);
}

int main(int argc, char* argv[]){
    // Get the number of processors
    int np = omp_get_num_procs();
    // Set the number of OpenMP threads to be the number of processors    
    omp_set_num_threads(np);

    
    // First, run it with the OpenMP turned on.
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double timing1 = tim.tv_sec+(tim.tv_usec/1000000.0);
    run_sim();

    // Reduce the number of threads to 1 and run again.
    gettimeofday(&tim, NULL);
    double timing2 = tim.tv_sec+(tim.tv_usec/1000000.0);
    omp_set_num_threads(1);
    run_sim();
    gettimeofday(&tim, NULL);
    double timing3 = tim.tv_sec+(tim.tv_usec/1000000.0);

    // Output speedup
    printf("\n\nOpenMP speed-up: %.3fx (perfect scaling would give %dx)\n",(timing3-timing2)/(timing2-timing1),np);
}

