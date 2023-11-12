/**
 * Shearing sheet with diagnostics
 *
 * This example simulates a small patch of Saturn's
 * Rings in shearing sheet coordinates. It also calculated
 * various quantities which can be used as diagnostics for
 * dynamical models of the rings. Diagnostics include
 * the midplane filling factor, the mean normal optical 
 * depth, the velocity dispersion tensor, the 
 * translational viscosity and the collisional viscosity.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

// This example is using a custom velocity dependend coefficient of restitution
double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v){
    // assumes v in units of [m/s]
    double eps = 0.32*pow(fabs(v)*100.,-0.234);
    if (eps>1) eps=1;
    if (eps<0) eps=0;
    return eps;
}

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->opening_angle2     = .5;                 // This determines the precission of the tree code gravity calculation.
    r->integrator         = REB_INTEGRATOR_SEI;
    r->boundary           = REB_BOUNDARY_SHEAR;
    r->gravity            = REB_GRAVITY_TREE;
    r->collision          = REB_COLLISION_TREE;
    r->collision_resolve  = reb_collision_resolve_hardsphere;
    double OMEGA          = 0.00013143527;      // 1/s
    r->ri_sei.OMEGA       = OMEGA;
    r->G                  = 6.67428e-11;        // N / (1e-5 kg)^2 m^2
    r->softening          = 0.1;                // m
    r->dt                 = 1e-3*2.*M_PI/OMEGA; // s
    r->heartbeat          = heartbeat;          // function pointer for callbacks after every timestep
    // This example uses two root boxes in the x and y direction. 
    // Although not necessary in this case, it allows for the parallelization using MPI. 
    // See Rein & Liu for a description of what a root box is in this context.
    double surfacedensity         = 400;        // kg/m^2
    double particle_density       = 400;        // kg/m^3
    double particle_radius_min    = 1;          // m
    double particle_radius_max    = 4;          // m
    double particle_radius_slope  = -3;    
    double boxsize                = 100;        // m
    reb_simulation_configure_box(r, boxsize, 2, 2, 1);
    r->N_ghost_x = 2;
    r->N_ghost_y = 2;
    r->N_ghost_z = 0;
    
    // Use Bridges et al coefficient of restitution.
    r->coefficient_of_restitution = coefficient_of_restitution_bridges;
    // When two particles collide and the relative velocity is zero, the might sink into each other in the next time step.
    // By adding a small repulsive velocity to each collision, we prevent this from happening.
    r->minimum_collision_velocity = particle_radius_min*OMEGA*0.001;  // small fraction of the shear accross a particle


    // Add all ring paricles
    double total_mass = surfacedensity*r->boxsize.x*r->boxsize.y;
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
    reb_simulation_integrate(r, INFINITY);
}

double mean_normal_geometric_optical_depth(const struct reb_simulation* const r){
    double area = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        area += M_PI*p.r*p.r;
    }
    return area/(r->boxsize.x*r->boxsize.y);
}

double midplane_fillingfactor(const struct reb_simulation* const r){
    double area = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double R2 = p.r*p.r-p.z*p.z;
        if (R2>0.){
            area += M_PI*R2;
        }
    }
    return area/(r->boxsize.x*r->boxsize.y);
}

struct reb_vec3d velocity_dispersion(const struct reb_simulation* const r){
    // Algorithm with reduced roundoff errors (see wikipedia)
    // Note: Average velocities relative to shear (stored in A) are not returned
    struct reb_vec3d A = {.x=0, .y=0, .z=0};
    struct reb_vec3d Q = {.x=0, .y=0, .z=0};
    for (int i=0;i<r->N;i++){
        struct reb_vec3d Aim1 = A;
        struct reb_particle p = r->particles[i];
        A.x = A.x + (p.vx-A.x)/(double)(i+1);
        A.y = A.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y)/(double)(i+1);
        A.z = A.z + (p.vz-A.z)/(double)(i+1);
        Q.x = Q.x + (p.vx-Aim1.x)*(p.vx-A.x);
        Q.y = Q.y + (p.vy+1.5*r->ri_sei.OMEGA*p.x-Aim1.y)*(p.vy+1.5*r->ri_sei.OMEGA*p.x-A.y);
        Q.z = Q.z + (p.vz-Aim1.z)*(p.vz-A.z);
    }
    Q.x = sqrt(Q.x/(double)r->N);
    Q.y = sqrt(Q.y/(double)r->N);
    Q.z = sqrt(Q.z/(double)r->N);

    // Return velocity dispersion in xx, yy, zz
    return Q; 
}

double translational_viscosity(const struct reb_simulation* const r){
    double Wxy = 0.;
    for (int i=0;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        double vx = p.vx;
        double vy = p.vy+1.5*r->ri_sei.OMEGA*p.x;
        Wxy += vx*vy;
    }
    return 2./3.*Wxy/r->N/r->ri_sei.OMEGA;
}

double collisional_viscosity(const struct reb_simulation* const r){
    // This is a time average!
    // To reset, set r->collisions_plog equal to 0.
    double Mtotal = 0.;
    for (int i=0;i<r->N;i++){
        Mtotal += r->particles[i].m;
    }
    return 2./3./r->ri_sei.OMEGA/Mtotal/r->t* r->collisions_plog;
}

void heartbeat(struct reb_simulation* const r){
    if (reb_simulation_output_check(r, 1e-3*2.*M_PI/r->ri_sei.OMEGA)){
        printf("Midplane FF=  %5.3f\t",midplane_fillingfactor(r));
        printf("Mean normal tau=  %5.3f \t",mean_normal_geometric_optical_depth(r));
        struct reb_vec3d Q = velocity_dispersion(r);
        printf("<vxvx>,<vyvy>,<vzvz>= %5.3e %5.3e %5.3e\t",Q.x, Q.y, Q.z);
        printf("nu_trans= %5.3e\t",translational_viscosity(r));
        printf("nu_col= %5.3e\t",collisional_viscosity(r));

        printf("\n");
    }
}

