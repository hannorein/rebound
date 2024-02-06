/**
 * Animation of the Saturn's Rings
 * 
 * This examples show how to use display_settings to
 * programmatically change the visualization of a 
 * REBOUND simulation. Here, we visualize a simulation of
 * Saturn's rings and rotate the viewwing angle programatically.
 * To understand what a 4x4 view matrix is, you can read 
 * up on linear algebra for computer graphics, specifically
 * the Model-View-Projection (MVP) paradigm. 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double coefficient_of_restitution_bridges(const struct reb_simulation* const r, double v);

// The heartbeat function is called once a timestep
void heartbeat(struct reb_simulation* const r){
    // Construct a rotation. 
    struct reb_vec3d axis = {.y=1., .z=0.2}; 
    struct reb_rotation rot = reb_rotation_init_angle_axis(0.003, axis); // small increment every timestep

    // Convert quaternion to rotation matrix
    struct reb_mat4df rm = reb_rotation_to_mat4df(rot);

    // Apply incremental rotation to view matrix.
    r->display_settings->view = reb_mat4df_multiply(rm, r->display_settings->view); 
}


int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();
    
    // The heartbeat function handles the visualization in this example.
    r->heartbeat         = heartbeat;           
    
    // Setup problem. For more details, see the shearing sheet example.
    r->opening_angle2    = .5;
    r->integrator        = REB_INTEGRATOR_SEI;
    r->boundary          = REB_BOUNDARY_SHEAR;
    r->gravity           = REB_GRAVITY_TREE;
    r->collision         = REB_COLLISION_TREE;
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->coefficient_of_restitution = coefficient_of_restitution_bridges;
    r->ri_sei.OMEGA      = 0.00013143527;       // 1/s
    r->minimum_collision_velocity = r->ri_sei.OMEGA*0.001;
    r->G                 = 6.67428e-11;         // N / (1e-5 kg)^2 m^2
    r->softening         = 0.1;                 // m
    r->dt                = 1e-3*2.*M_PI/r->ri_sei.OMEGA;  // s
    
    reb_simulation_configure_box(r, 100, 2, 2, 1); // 100m box
    r->N_ghost_x = 2; r->N_ghost_y = 2; r->N_ghost_z = 0;
    
    // Add all ring paricles
    double mass = 0;
    while(mass < 400*r->boxsize.x*r->boxsize.y){ // 400kg/m^2 surface density
        struct reb_particle pt = {0};
        pt.x         = reb_random_uniform(r, -r->boxsize.x/2.,r->boxsize.x/2.);
        pt.y         = reb_random_uniform(r, -r->boxsize.y/2.,r->boxsize.y/2.);
        pt.z         = reb_random_normal(r, 1.);                    // m
        pt.vy         = -1.5*pt.x*r->ri_sei.OMEGA;
        double radius     = reb_random_powerlaw(r, 1., 4.,-3);
        pt.r         = radius;                        // m
        double        particle_mass = 400.0*4./3.*M_PI*radius*radius*radius; // 400kg/m^3 particle density
        pt.m         = particle_mass;     // kg
        reb_simulation_add(r, pt);
        mass += particle_mass;
    }
    
    // Normally the visualization settings are determined by the 
    // user interface. If we add the display_settings struct to 
    // the simulation itself, it will overwrite any change the 
    // user has made and allows us to programatically change any 
    // settings such as the orientation, zoom, etc. 
    reb_simulation_add_display_settings(r);

    // This allows you to connect to the simulation using
    // a web browser. Simply go to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Integrate forever
    reb_simulation_integrate(r, INFINITY);
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

