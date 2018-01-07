/**
 * Granular dynamics.
 *
 * This example is about granular dynamics. No gravitational 
 * forces are present in this example. Two boundary layers made of 
 * particles simulate shearing walls. These walls are heating
 * up the particles, create a dense and cool layer in the middle.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

int collision_resolve_hardsphere_withborder(struct reb_simulation* r, struct reb_collision c);
void heartbeat(struct reb_simulation* r);
int N_border;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->dt             = 1e-1;    
    r->gravity        = REB_GRAVITY_NONE;
    r->integrator        = REB_INTEGRATOR_LEAPFROG;
    r->collision        = REB_COLLISION_TREE;
    r->boundary        = REB_BOUNDARY_PERIODIC;
    // Override default collision handling to account for border particles
    r->collision_resolve     = collision_resolve_hardsphere_withborder;
    r->heartbeat        = heartbeat;
    reb_configure_box(r, 20., 1, 1, 4);
    
    r->nghostx = 1; r->nghosty = 1; r->nghostz = 0;     
    
    double N_part     = 0.00937*r->boxsize.x*r->boxsize.y*r->boxsize.z;

    // Add Border Particles
    double radius         = 1;
    double mass        = 1;
    double border_spacing_x = r->boxsize.x/(floor(r->boxsize.x/radius/2.)-1.);
    double border_spacing_y = r->boxsize.y/(floor(r->boxsize.y/radius/2.)-1.);
    struct reb_particle pt = {0};
    pt.r         = radius;
    pt.m         = mass;
    pt.hash        = 1;
    for(double x = -r->boxsize.x/2.; x<r->boxsize.x/2.-border_spacing_x/2.;x+=border_spacing_x){
        for(double y = -r->boxsize.y/2.; y<r->boxsize.y/2.-border_spacing_y/2.;y+=border_spacing_y){
            pt.x         = x;
            pt.y         = y;
            
            // Add particle to bottom
            pt.z         = -r->boxsize.z/2.+radius;
            pt.vy         = 1;
            reb_add(r, pt);

            // Add particle to top
            pt.z         = r->boxsize.z/2.-radius;
            pt.vy         = -1;
            reb_add(r, pt);
        }
    }

    N_border = r->N;
    
    // Add real particles
    while(r->N-N_border<N_part){
        struct reb_particle pt = {0};
        pt.x         = reb_random_uniform(-r->boxsize.x/2.,r->boxsize.x/2.);
        pt.y         = reb_random_uniform(-r->boxsize.y/2.,r->boxsize.y/2.);
        pt.z         = 0.758*reb_random_uniform(-r->boxsize.z/2.,r->boxsize.z/2.);
        pt.vx         = reb_random_normal(0.001);
        pt.vy         = reb_random_normal(0.001);
        pt.vz         = reb_random_normal(0.001);
        pt.r         = radius;                        // m
        pt.m         = 1;
        pt.hash        = 2;
        reb_add(r, pt);
    }

    reb_integrate(r, INFINITY);
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10.*r->dt)){
        reb_output_timing(r, 0);
    }
}

int collision_resolve_hardsphere_withborder(struct reb_simulation* r, struct reb_collision c){
    const double t = r->t;
    struct reb_particle* particles = r->particles;
    struct reb_particle p1 = particles[c.p1];
    struct reb_particle p2 = particles[c.p2];
    struct reb_ghostbox gb = c.gb;
    double m21  = p1.m  /  p2.m; 
    double x21  = p1.x + gb.shiftx  - p2.x; 
    double y21  = p1.y + gb.shifty  - p2.y; 
    double z21  = p1.z + gb.shiftz  - p2.z; 
    double rp   = p1.r+p2.r;
    if (rp*rp < x21*x21 + y21*y21 + z21*z21) return 0;
    double vx21 = p1.vx + gb.shiftvx - p2.vx; 
    double vy21 = p1.vy + gb.shiftvy - p2.vy; 
    double vz21 = p1.vz + gb.shiftvz - p2.vz; 
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

    // Coefficient of restitution
    double eps = 0.15;
    double dvx2 = -(1.0+eps)*vx21nn/(1.0+m21) ;

    // Now we are rotating backwards
    double dvx2n = cphi * dvx2;        
    double dvy2n = sphi * dvx2;        
    double dvy2nn = ctheta * dvy2n;    
    double dvz2nn = stheta * dvy2n;    

    // Applying the changes to the particles.
    // Do not change border particles.
    if (p2.hash!=1){
        particles[c.p2].vx -=    m21*dvx2n;
        particles[c.p2].vy -=    m21*dvy2nn;
        particles[c.p2].vz -=    m21*dvz2nn;
        particles[c.p2].lastcollision = t;
    }
    if (p1.hash!=1){
        particles[c.p1].vx +=    dvx2n; 
        particles[c.p1].vy +=    dvy2nn; 
        particles[c.p1].vz +=    dvz2nn; 
        particles[c.p1].lastcollision = t;
    }
    return 0; // Do not remove any particle from simulation.
}
