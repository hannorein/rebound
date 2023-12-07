/**
 * Shearing sheet with MPI
 *
 * This example simulates a small patch of Saturn's
 * Rings in shearing sheet coordinates. The code can use MPI 
 * to distribute the work of the gravity and collision modules
 * to other nodes. You can enable OpenGL with MPI, but this 
 * is a feature that might not work in all environments. 
 * You can turn on OpenGL in the Makefile.
 * How to configure and submit an MPI job varies significantly
 * depending on your cluster architecture. To test MPI on your
 * local computer, simply type make && mpirun -np 4 rebound.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v);
void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->opening_angle2    = .5;                    // This determines the precission of the tree code gravity calculation.
    r->integrator        = REB_INTEGRATOR_SEI;
    r->boundary          = REB_BOUNDARY_SHEAR;
    r->gravity           = REB_GRAVITY_TREE;
    r->collision         = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    double OMEGA         = 0.00013143527;    // 1/s
    r->ri_sei.OMEGA      = OMEGA;
    r->G                 = 6.67428e-11;        // N / (1e-5 kg)^2 m^2
    r->softening         = 0.1;            // m
    r->dt                = 1e-3*2.*M_PI/OMEGA;    // s
    r->heartbeat         = heartbeat;    // function pointer for heartbeat
    // This example uses two root boxes in the x and y direction. 
    // Although not necessary in this case, it allows for the parallelization using MPI. 
    // See Rein & Liu for a description of what a root box is in this context.
    double surfacedensity          = 400;          // kg/m^2
    double particle_density        = 400;          // kg/m^3
    double particle_radius_min     = 1;            // m
    double particle_radius_max     = 4;            // m
    double particle_radius_slope     = -3;    
    double boxsize             = 100;              // m
    if (argc>1){                        // Try to read boxsize from command line
        boxsize = atof(argv[1]);
    }
    // Setup 2x2 root boxes.
    // This allows you to use up to 4 MPI nodes.
    reb_simulation_configure_box(r, boxsize, 2, 2, 1);
    // Initialize MPI (this only works after reb_simulation_configure_box)
    reb_mpi_init(r);
    r->N_ghost_x = 2;
    r->N_ghost_y = 2;
    r->N_ghost_z = 0;
    
    // Initial conditions
    printf("Toomre wavelength: %f\n",4.*M_PI*M_PI*surfacedensity/OMEGA/OMEGA*r->G);
    // Use Bridges et al coefficient of restitution.
    r->coefficient_of_restitution = coefficient_of_restitution_bridges;
    // When two particles collide and the relative velocity is zero, the might sink into each other in the next time step.
    // By adding a small repulsive velocity to each collision, we prevent this from happening.
    r->minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear accross a particle


    // Add all ring paricles
    double total_mass = surfacedensity*r->boxsize.x*r->boxsize.y/r->mpi_num;
    double mass = 0;
    while(mass<total_mass){
        struct reb_particle pt;
        pt.x         = reb_random_uniform(r, -r->boxsize.x/2.,r->boxsize.x/2.);
        pt.y         = reb_random_uniform(r, -r->boxsize.y/2.,r->boxsize.y/2.);
        pt.z         = reb_random_normal(r, 1.);                    // m
        pt.vx         = 0;
        pt.vy         = -1.5*pt.x*OMEGA;
        pt.vz         = 0;
        pt.ax         = 0;
        pt.ay         = 0;
        pt.az         = 0;
        double radius     = reb_random_powerlaw(r, particle_radius_min,particle_radius_max,particle_radius_slope);
        pt.r         = radius;                        // m
        double        particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
        pt.m         = particle_mass;     // kg
        reb_simulation_add(r, pt);
        mass += particle_mass;
    }
#ifdef OPENGL
    // Hack to artificially increase particle array.
    // This cannot be done once OpenGL is activated. 
    r->N_allocated *=8;
    r->particles = realloc(r->particles,sizeof(struct reb_particle)*r->N_allocated);
#endif // OPENGL

    // Start the integration
    reb_simulation_integrate(r, INFINITY);

    // Cleanup
    reb_mpi_finalize(r);
    reb_simulation_free(r);
}

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
    // assumes v in units of [m/s]
    double eps = 0.32*pow(fabs(v)*100.,-0.234);
    if (eps>1) eps=1;
    if (eps<0) eps=0;
    return eps;
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 1e-3*2.*M_PI/r->ri_sei.OMEGA)){
        reb_simulation_output_timing(r, 0);
        //reb_output_append_velocity_dispersion("veldisp.txt");
    }
    if (reb_simulation_output_check(r, 2.*M_PI/r->ri_sei.OMEGA)){
        //reb_simulation_output_ascii("position.txt");
    }
}

