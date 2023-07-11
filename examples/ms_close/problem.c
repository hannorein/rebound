/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
    struct reb_particle star = {0};
    star.m = 1;

    // Jupiter
    double jm = 9.55e-4;
    double ja = 5.2;
    double je = 0.05;

    // Saturn
    double sm = 2.857e-4;
    double sa = 9.58;
    double se = 0.95;
    double si = M_PI / 2.;

    r->dt = 0.1*2.*M_PI;
    r->integrator = REB_INTEGRATOR_BS;

    //r->integrator = REB_INTEGRATOR_TRACE;
    //r->ri_tr.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    //r->ri_tr.peri = 100;
    r->heartbeat  = heartbeat;

    reb_add(r, star);
    reb_add_fmt(r, "m a e", jm, ja, je);

    //struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
    //r->dt = orb.P / 15.12345;
    reb_add_fmt(r, "m a e inc omega", sm, sa, se, si, M_PI/2);

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_tools_energy(r);
    system("rm -rf energy.txt");
    //FILE* f = fopen("energy.txt","w");

    reb_integrate(r, 1000*2*M_PI*11.86);
    //reb_integrate(r, 1277.);
    //err = reb_tools_energy(r);
    //printf("%e\n", (err - e_init) /  e_init);
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}
    if (reb_output_check(r, 1. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");

        // rotate whole simulation to rotating frame
        //reb_simulation_irotate(r, r1);
        struct reb_particle* sun = &r->particles[0];
        struct reb_particle* jup = &r->particles[1];
        struct reb_particle* sat = &r->particles[2];

        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *jup, *sun);
        double e = reb_tools_energy(r);
        fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %f %f\n",r->t, (e - e_init) / e_init, r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z, o.e, o.inc);
        fclose(f);

        //reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);
    }
}
