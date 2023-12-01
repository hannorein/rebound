/**
 * Planetary migration in the GJ876 system
 *
 * This example applies dissipative forces to two
 * bodies orbiting a central object. The forces are specified
 * in terms of damping timescales for the semi-major axis and
 * eccentricity. This mimics planetary migration in a protostellar disc. 
 * The example reproduces the study of Lee & Peale (2002) on the 
 * formation of the planetary system GJ876. For a comparison, 
 * see figure 4 in their paper. The IAS15 or WHFAST integrators
 * can be used. Note that the forces are velocity dependent.
 * Special thanks goes to Willy Kley for helping me to implement
 * the damping terms as actual forces. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

double* tau_a;     /**< Migration timescale in years for all particles */
double* tau_e;     /**< Eccentricity damping timescale in years for all particles */
double tmax;

void migration_forces(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    // Start the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);

    // Setup constants
    r->integrator           = REB_INTEGRATOR_WHFAST;
    //r->integrator         = REB_INTEGRATOR_IAS15;
    r->dt                   = 1e-2*2.*M_PI;        // in year/(2*pi)
    r->additional_forces    = migration_forces;     //Set function pointer to add dissipative forces.
    r->heartbeat            = heartbeat;  
    r->force_is_velocity_dependent = 1;
    tmax                    = 2.0e4*2.*M_PI;    // in year/(2*pi)

    // Initial conditions
    // Parameters are those of Lee & Peale 2002, Figure 4. 
    struct reb_particle star = {0};
    star.m  = 0.32;            // This is a sub-solar mass star
    reb_simulation_add(r, star); 
    
    struct reb_particle p1 = {0};    // Planet 1
    p1.x     = 0.5;
    p1.m      = 0.56e-3;
    p1.vy     = sqrt(r->G*(star.m+p1.m)/p1.x);
    reb_simulation_add(r, p1); 
    
    struct reb_particle p2 = {0};    // Planet 2
    p2.x     = 1;
    p2.m      = 1.89e-3;
    p2.vy     = sqrt(r->G*(star.m+p2.m)/p2.x);
    reb_simulation_add(r, p2); 

    tau_a = calloc(sizeof(double),r->N);
    tau_e = calloc(sizeof(double),r->N);

    tau_a[2] = 2.*M_PI*20000.0;    // Migration timescale of planet 2 is 20000 years.
    tau_e[2] = 2.*M_PI*200.0;      // Eccentricity damping timescale is 200 years (K=100). 

    reb_simulation_move_to_com(r);          

    remove("orbits.txt"); // delete previous output file

    reb_simulation_integrate(r, tmax);
}

void migration_forces(struct reb_simulation* r){
    const double G = r->G;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0]; // calculate migration forces with respect to center of mass;
    for(int i=1;i<N;i++){
        if (tau_e[i]!=0||tau_a[i]!=0){
            struct reb_particle* p = &(particles[i]);
            const double dvx = p->vx-com.vx;
            const double dvy = p->vy-com.vy;
            const double dvz = p->vz-com.vz;

            if (tau_a[i]!=0){     // Migration
                p->ax -=  dvx/(2.*tau_a[i]);
                p->ay -=  dvy/(2.*tau_a[i]);
                p->az -=  dvz/(2.*tau_a[i]);
            }
            if (tau_e[i]!=0){     // Eccentricity damping
                const double mu = G*(com.m + p->m);
                const double dx = p->x-com.x;
                const double dy = p->y-com.y;
                const double dz = p->z-com.z;

                const double hx = dy*dvz - dz*dvy; 
                const double hy = dz*dvx - dx*dvz;
                const double hz = dx*dvy - dy*dvx;
                const double h = sqrt ( hx*hx + hy*hy + hz*hz );
                const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
                const double r = sqrt ( dx*dx + dy*dy + dz*dz );
                const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
                const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
                const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
                const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
                const double e = sqrt( ex*ex + ey*ey + ez*ez );        // eccentricity
                const double a = -mu/( v*v - 2.*mu/r );            // semi major axis
                const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
                const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
                p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
                p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
                p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
            }
        }
        com = reb_particle_com_of_pair(com,particles[i]);
    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r, 20.*M_PI)){
        reb_simulation_output_timing(r, tmax);
    }
    if(reb_simulation_output_check(r, 40.)){
        reb_simulation_synchronize(r);
        reb_simulation_output_orbits(r,"orbits.txt");
        reb_simulation_move_to_com(r); 
    }
}
