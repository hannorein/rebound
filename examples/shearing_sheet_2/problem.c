/**
 * Shearing sheet (Akihiko Fujii)
 *
 * This example is identical to the shearing_sheet
 * example but uses a different algorithm for resolving individual 
 * collisions. In some cases, this might give more realistic results.
 * Particle properties resemble those found in Saturn's rings. 
 *
 * In this collision resolve method, particles are displaced if they 
 * overlap. This example also shows how to implement your own collision
 * routine. This is where one could add fragmentation, or merging of
 * particles.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

int collision_resolve_hardsphere_pullaway(struct reb_simulation* r, struct reb_collision c);

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v);
void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->opening_angle2    = .5;                    // This determines the precission of the tree code gravity calculation.
    r->integrator        = REB_INTEGRATOR_SEI;
    r->boundary          = REB_BOUNDARY_SHEAR;
    r->gravity           = REB_GRAVITY_TREE;
    r->collision         = REB_COLLISION_TREE;
    r->collision_resolve = collision_resolve_hardsphere_pullaway;
    double OMEGA         = 0.00013143527;       // 1/s
    r->ri_sei.OMEGA      = OMEGA;
    r->G                 = 6.67428e-11;         // N / (1e-5 kg)^2 m^2
    r->softening         = 0.1;                 // m
    r->dt                = 1e-3*2.*M_PI/OMEGA;  // s
    r->heartbeat         = heartbeat;           // function pointer for heartbeat
    // This example uses two root boxes in the x and y direction. 
    // Although not necessary in this case, it allows for the parallelization using MPI. 
    // See Rein & Liu for a description of what a root box is in this context.
    double surfacedensity          = 400;      // kg/m^2
    double particle_density        = 400;      // kg/m^3
    double particle_radius_min     = 1;        // m
    double particle_radius_max     = 4;        // m
    double particle_radius_slope   = -3;    
    double boxsize                 = 100;      // m
    if (argc>1){                               // Try to read boxsize from command line
        boxsize = atof(argv[1]);
    }
    reb_simulation_configure_box(r, boxsize, 2, 2, 1);
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
    double total_mass = surfacedensity*r->boxsize.x*r->boxsize.y;
    double mass = 0;
    while(mass<total_mass){
        struct reb_particle pt = {0};
        pt.x         = reb_random_uniform(r, -r->boxsize.x/2.,r->boxsize.x/2.);
        pt.y         = reb_random_uniform(r, -r->boxsize.y/2.,r->boxsize.y/2.);
        pt.z         = reb_random_normal(r, 1.);                    // m
        pt.vy         = -1.5*pt.x*OMEGA;
        double radius     = reb_random_powerlaw(r, particle_radius_min,particle_radius_max,particle_radius_slope);
        pt.r         = radius;                        // m
        double        particle_mass = particle_density*4./3.*M_PI*radius*radius*radius;
        pt.m         = particle_mass;     // kg
        reb_simulation_add(r, pt);
        mass += particle_mass;
    }
    reb_simulation_integrate(r, INFINITY);
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

// Function written by Akihiko Fujii
int collision_resolve_hardsphere_pullaway(struct reb_simulation* r, struct reb_collision c){
    struct reb_particle* particles = r->particles;
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2 = particles[c.p2];
    struct reb_vec6d gb = c.gb;
    double x21  = p1.x + gb.x  - p2.x; 
    double y21  = p1.y + gb.y  - p2.y; 
    double z21  = p1.z + gb.z  - p2.z; 
    double _r = sqrt(x21*x21 + y21*y21 + z21*z21);
    /* double r21 = sqrt(x21*x21 + y21*y21 + z21*z21); */
    double rp   = p1.r+p2.r;
    
    if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;

    double vx21 = p1.vx + gb.vx - p2.vx; 
    double vy21 = p1.vy + gb.vy - p2.vy; 
    double vz21 = p1.vz + gb.vz - p2.vz; 

    if (vx21*x21 + vy21*y21 + vz21*z21 >0) return 0; // not approaching

    // Bring the to balls in the xy plane.
    // NOTE: this could probabely be an atan (which is faster than atan2) 
    double theta = atan2(z21,y21);
    double stheta = sin(theta);
    double ctheta = cos(theta);
    double vy21n = ctheta * vy21 + stheta * vz21;
    double y21n = ctheta * y21 + stheta * z21;

    // Bring the two balls onto the positive x axis.
    double phi = atan2(y21n,x21);
    double cphi = cos(phi);
    double sphi = sin(phi);
    double vx21nn = cphi * vx21  + sphi * vy21n;
    double vy21nn = -sphi* vx21  + cphi * vy21n;

    // Coefficient of restitution
    double eps= r->coefficient_of_restitution(r, vx21nn);
    double dvx2 = -(1.0+eps)*vx21nn;
    double dvy2 = (_r/rp-1.)*vy21nn;

    double minr = (p1.r>p2.r)?p2.r:p1.r;
    double maxr = (p1.r<p2.r)?p2.r:p1.r;
    double mindv= minr*r->minimum_collision_velocity;
    mindv *= 1.-(_r - maxr)/minr;
    if (mindv>maxr*r->minimum_collision_velocity)mindv = maxr*r->minimum_collision_velocity;
    if (dvx2<mindv) dvx2 = mindv;

    double dxx2 = rp-_r;
    double dxx2n = cphi * dxx2;
    double dxy2n = sphi * dxx2;
    double dxy2nn = ctheta * dxy2n;
    double dxz2nn = stheta * dxy2n;

    double dvx2n = cphi * dvx2 - sphi * dvy2;
    double dvy2n = sphi * dvx2 + cphi * dvy2;

    double dvy2nn = ctheta * dvy2n;
    double dvz2nn = stheta * dvy2n;

    // Applying the changes to the particles.
    const double p1pf = p1.m/(p1.m+p2.m);
    const double p2pf = p2.m/(p1.m+p2.m);
    particles[c.p2].vx -=    p1pf*dvx2n;
    particles[c.p2].vy -=    p1pf*dvy2nn;
    particles[c.p2].vz -=    p1pf*dvz2nn;
    particles[c.p2].last_collision = r->t;
    particles[c.p2].x -=    p1pf*dxx2n;
    particles[c.p2].y -=    p1pf*dxy2nn;
    particles[c.p2].z -=    p1pf*dxz2nn;


    particles[c.p1].vx +=    p2pf*dvx2n;
    particles[c.p1].vy +=    p2pf*dvy2nn;
    particles[c.p1].vz +=    p2pf*dvz2nn;
    particles[c.p1].x +=    p2pf*dxx2n; 
    particles[c.p1].y +=    p2pf*dxy2nn; 
    particles[c.p1].z +=    p2pf*dxz2nn; 

    particles[c.p1].last_collision = r->t;

    return 0; // Do not remove any particle
}

