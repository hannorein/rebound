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
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->ri_mercurius.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    r->ri_mercurius.peri = 1;
    r->heartbeat  = heartbeat;

    reb_add(r, star);
    reb_add_fmt(r, "m a e", jm, ja, je);

    //struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[1], r->particles[0]);
    //r->dt = orb.P / 15.12345;
    reb_add_fmt(r, "m a e inc f", sm, sa, se, si, M_PI);

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_tools_energy(r);
    system("rm -rf energy.txt");
    FILE* f = fopen("energy.txt","w");

    reb_integrate(r, 3000 * 2 * M_PI);

    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}
    if (reb_output_check(r, (4. / 365.25) * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");

        // rotate whole simulation to rotating frame
        //reb_simulation_irotate(r, r1);
        double e = reb_tools_energy(r);
        fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e\n",r->t, (e - e_init) / e_init, r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z);
        fclose(f);

        //reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);
    }
}
