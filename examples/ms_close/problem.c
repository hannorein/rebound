/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;
double emax = 0.0;
double tmax = 300. * 2. * M_PI * 29.4;

const double k = 1.; // Constants for the Harmonic Oscillator
const double m = 10000000.;

void derivatives(struct reb_ode* const ode, double* const yDot, const double* const y, const double t){
    const double omega = sqrt(k/m);
    //struct reb_orbit o = reb_orbit_from_particle(ode->r->G, ode->r->particles[1], ode->r->particles[0]);
    //double forcing = sin(o.f);
    yDot[0] = y[1];
    yDot[1] = -omega*omega*y[0];// + forcing;
}

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle star = {0};
    star.m = 1;

    // Jupiter
    double jm = 9.55e-4;
    double ja = 5.2;
    double je = 0.05;
    double jmse = 1.-je;

    // Saturn
    double sm = 2.857e-4;
    double sa = 9.58;
    double se = 0.99;
    double omse = 1. - se;
    double si = M_PI / 2.;

    //const double tau_f_j = 2 * M_PI * sqrt((ja * ja * ja * jmse * jmse * jmse) / (r->G * star.m * (1 + je)));
    //const double tau_f_s = 2 * M_PI * sqrt((sa * sa * sa * omse * omse * omse) / (r->G * star.m * (1 + se)));

    r->integrator = REB_INTEGRATOR_TRACE;
    r->dt = 0.15 * 2 * M_PI;
    //r->ri_bs.eps_rel = 1e-11;            // Relative tolerance
    //r->ri_bs.eps_abs = 1e-11;            // Absolute tolerance
    //r->ri_ias15.adaptive_mode=2;
    r->ri_trace.peri_crit_fdot=16.;

    //struct reb_ode* ho = reb_ode_create(r,2);   // Add an ODE with 2 dimensions
    //ho->derivatives = derivatives;              // Right hand side of the ODE
    //ho->y[0] = 1;                               // Initial conditions
    //ho->y[1] = 0;

    //r->dt = 0.15*2.*M_PI;
    r->heartbeat = heartbeat;


    reb_simulation_add(r, star);
    reb_simulation_add_fmt(r, "m a e", jm, ja, je);

    //r->dt = orb.P / 15.12345;
    reb_simulation_add_fmt(r, "primary m a e inc omega", star, sm, sa, se, si, M_PI/2.);


    struct reb_orbit orb = reb_orbit_from_particle(r->G, r->particles[2], r->particles[0]);

    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    //r->heartbeat  = heartbeat;
    r->exact_finish_time=0;
    e_init = reb_simulation_energy(r);
    system("rm -rf ho.txt");
    system("rm -rf energy_fdot.txt");

    FILE* f = fopen("ho.txt","a");
    while(r->t<tmax){
        reb_simulation_integrate(r, r->t + (tmax/3000.));
        //fprintf(f, "%f,%f\n", r->t, ho->y[0]);
    }
    fclose(f);
    printf("\n%e\n", fabs((reb_simulation_energy(r) - e_init) / e_init));
    reb_simulation_free(r);
}


void heartbeat(struct reb_simulation* r){
  //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
    if (reb_simulation_output_check(r, 10.*2.*M_PI)){
        reb_simulation_output_timing(r, tmax);
    }

    //if (reb_output_check(r, 10. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file

        double err = fabs((reb_simulation_energy(r) - e_init) / e_init);
        if (err > emax){
          emax = err;
        }
        FILE* f = fopen("energy_fdot.txt","a");

        // rotate whole simulation to rotating frame
        //reb_simulation_irotate(r, r1);

        struct reb_particle* sun = &r->particles[0];
        struct reb_particle* jup = &r->particles[1];
        struct reb_particle* sat = &r->particles[2];

        struct reb_orbit o = reb_orbit_from_particle(r->G, *sat, *sun);
        struct reb_orbit oj = reb_orbit_from_particle(r->G, *jup, *sun);

        const double dvx = r->particles[0].vx - r->particles[2].vx;
        const double dvy = r->particles[0].vy - r->particles[2].vy;
        const double dvz = r->particles[0].vz - r->particles[2].vz;
        const double dv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);

        const double dvxj = r->particles[0].vx - r->particles[1].vx;
        const double dvyj = r->particles[0].vy - r->particles[1].vy;
        const double dvzj = r->particles[0].vz - r->particles[1].vz;
        const double dvj = sqrt(dvxj*dvxj + dvyj*dvyj + dvzj*dvzj);

        const double vcirc = sqrt(r->G * r->particles[0].m / o.a);
        const double vcircj = sqrt(r->G * r->particles[0].m / oj.a);
        double e = (reb_simulation_energy(r) - e_init) / e_init;
        fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %f %f %f %f\n",r->t, e, r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z, o.e, o.inc, dvx/vcirc, dvxj/vcircj);

        //double e = reb_tools_energy(r);
        //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
        fclose(f);

        //reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);
    //}


}
