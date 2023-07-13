/**
 * Planetesimal Disk Migration
 *
 * This example integrates a star, 2 planet, N-planetesimal disk system, with the
 * outer planet at the inner edge of the planetesimal disk. The planet in the system
 * migrates on a very slow timescale and one needs to run the simulation for roughly
 * 10^5 dynamical timescales to see the effect.
 *
 * The ideal integrator choice for this problem is MERCURIUS due to the large
 * number of close encounters.  
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);
double E0;

int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c);

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    // Simulation Setup
    r->integrator    = REB_INTEGRATOR_MERCURIUS;
    r->heartbeat    = heartbeat;
    // Test particle type 1 allows massive particles to feel the gravity of testparticles.
    // However, test particles will not feel the gravity from other test particles.
    r->testparticle_type = 1;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge_pass_through;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary    = REB_BOUNDARY_OPEN;
    const double boxsize = 6;
    reb_configure_box(r,boxsize,2,2,1);
    
    srand(12);
    double m_earth = 3.003e-6;
    double m_neptune = 5.1e-4;
    double a_scat_planet = 1;
    double a_mig_planet = 1.67;
    r->dt = 6.283*pow(a_scat_planet,1.5)/50;
    
    // Star
    struct reb_particle star = {0};
    star.m         = 1;
    star.r        = 0.005;        // Radius of particle is in AU!
    reb_add(r, star);
    
    // Planet 1 - inner massive planet to scatter planetesimals out
    {
        double a=a_scat_planet, m=m_neptune, e=0, inc=reb_random_normal(r, 0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p.r = 0.000467;
        reb_add(r, p);
    }
    
    // Planet 2 - outer smaller planet to migrate in the disk
    {
        double a=a_mig_planet, m=2.3*m_earth, e=0, inc=reb_random_normal(r, 0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, e, inc, 0, 0, 0);
        p.r = 0.0000788215;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    
    // Planetesimal disk parameters
    double total_disk_mass = 2.3*10*m_earth;
    int N_planetesimals = 2500;
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    double amin = a_mig_planet-0.02, amax = a_mig_planet + 1;   //planet at inner edge of disk
    double powerlaw = 1;
    
    // Generate Planetesimal Disk
    while(r->N<N_planetesimals + r->N_active){
        struct reb_particle pt = {0};
        double a    = reb_random_powerlaw(r, amin,amax,powerlaw);
        double e    = reb_random_rayleigh(r, 0.005);
        double inc  = reb_random_rayleigh(r, 0.005);
        double Omega = reb_random_uniform(r, 0,2.*M_PI);
        double apsis = reb_random_uniform(r, 0,2.*M_PI);
        double phi     = reb_random_uniform(r, 0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
        pt.r         = 0.00000934532;
        reb_add(r, pt);
    }

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
    // Integrate!
    reb_integrate(r, INFINITY);
    reb_free_simulation(r);
}

int reb_collision_resolve_merge_pass_through(struct reb_simulation* const r, struct reb_collision c){
    // This function passes the collision to the default merging routine. 
    // If a merger occured, that routine will return a value other than 0.
    // This function then outputs some information about the merger.
    int result = reb_collision_resolve_merge(r, c);
    if (result!=0){
        printf("A merger occured! Particles involved: %d, %d.\n",c.p1,c.p2);
    }
    return result;
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 100.*r->dt)){
        //relative energy error
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        
        //get orbital elements
        struct reb_particle p = r->particles[2];
        struct reb_particle star = r->particles[0];
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G,p,star);
        
        printf("a2=%f,dE=%e,N=%d\n",o.a,relE,r->N);
    }
}
