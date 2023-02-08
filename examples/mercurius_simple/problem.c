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

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    r->dt = 0.1*2.*M_PI;
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->ri_mercurius.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
    r->heartbeat  = heartbeat;

    // Initialize masses
    struct reb_particle star = {0};
    star.m = 1;
    struct reb_particle jup = {0};
    jup.m = 0.01 / (star.m - 0.01);

    // velocities
    double a = 5.2;
    double e = 0.0;
    star.x = -(jup.m / (star.m + jup.m)) * (a * (1 + e));
    star.vy = -(jup.m / (star.m + jup.m)) * sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
    reb_add(r, star);

    jup.x = (star.m / (star.m + jup.m)) * (a * (1 + e));
    jup.vy = (star.m / (star.m + jup.m)) * sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
    reb_add(r, jup);

    // Test particle
    struct reb_particle test = {0};
    double xhel = 4.42;
    double vhel = 0.0072;

    test.x = xhel + star.x;
    test.vy = vhel + star.vy;
    reb_add(r, test);


    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_tools_energy(r);
    system("rm -rf energy.txt");

    reb_integrate(r, 50.*12.*2.*M_PI);
}

void heartbeat(struct reb_simulation* r){
    if (reb_output_check(r, 10.*2.*M_PI)){
        reb_output_timing(r, 0);
    }
    if (reb_output_check(r, 2.*M_PI)){
        // Once per year, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");
        reb_integrator_synchronize(r);
        double e = reb_tools_energy(r);
        fprintf(f,"%e %e\n",r->t, fabs((e-e_init)/e_init));
        fclose(f);
    }
}
