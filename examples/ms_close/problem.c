#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include <time.h>

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;
double emax = 0.0;
double tmax = 3000. * 2. * M_PI;

char title_stats[100] = "peri_stats";

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

    double peri = 0.1;
    int integrator = 0;
    if (argc == 3){
      peri = pow(10, -1. * atof(argv[2]));
      integrator = atoi(argv[1]);
    }

    // Saturn
    double sm = 2.857e-4;
    double sa = 9.58;
    double se = (sa - peri) / sa;
    double omse = 1. - se;
    double si = M_PI / 2.;

    char grators[10][10] = {"TRACE", "WHFAST", "MERCURIUS", "WHFASTr", "MERCURIUSr", "TRACEf"};

    switch(integrator){
      case 0:
        r->integrator = REB_INTEGRATOR_TRACE;
        r->dt = 0.15 * 2 * M_PI;
        r->ri_trace.peri_crit_distance=2.;
      case 1:
        r->integrator = REB_INTEGRATOR_WHFAST;
        r->dt = 0.15 * 2 * M_PI;
      case 2:
        r->integrator = REB_INTEGRATOR_MERCURIUS;
        r->dt = 0.15 * 2 * M_PI;
      case 3:
        r->integrator = REB_INTEGRATOR_WHFAST;
        r->dt = 2*M_PI*sqrt((1-se)*(1-se)*(1-se)*sa*sa*sa/((1+se)*r->G*star.m))/50.;
      case 4:
        r->integrator = REB_INTEGRATOR_MERCURIUS;
        r->dt = 2*M_PI*sqrt((1-se)*(1-se)*(1-se)*sa*sa*sa/((1+se)*r->G*star.m))/15.;
      case 5:
        r->integrator = REB_INTEGRATOR_TRACE;
        r->dt = 0.15 * 2 * M_PI;
        r->ri_trace.peri_crit_fdot=17.;
    }

    //r->integrator = REB_INTEGRATOR_TRACE;
    //r->dt = 2*M_PI*sqrt((1-se)*(1-se)*(1-se)*sa*sa*sa/((1+se)*r->G*star.m))/15.;//0.15 * 2 * M_PI;
    //r->dt = 0.15 * 2 * M_PI;
    //r->ri_trace.r_crit_hill = 0.;
    //r->ri_ias15.adaptive_mode=2;
    //r->ri_trace.peri_crit_fdot=0.;

    //struct reb_ode* ho = reb_ode_create(r,2);   // Add an ODE with 2 dimensions
    //ho->derivatives = derivatives;              // Right hand side of the ODE
    //ho->y[0] = 1;                               // Initial conditions
    //ho->y[1] = 0;
    reb_simulation_add(r, star);
    reb_simulation_add_fmt(r, "m a e", jm, ja, je);

    //r->dt = orb.P / 15.12345;
    reb_simulation_add_fmt(r, "primary m a e inc omega", star, sm, sa, se, si, M_PI/2.);


    struct reb_orbit orb = reb_orbit_from_particle(r->G, r->particles[2], r->particles[0]);

    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    r->heartbeat  = heartbeat;
    r->exact_finish_time=0;
    e_init = reb_simulation_energy(r);
    //system("rm -rf ho.txt");
    //system("rm -rf energy_trace_best.txt");

    //FILE* f = fopen("energy_trace_best.txt","a");
    //fclose(f);

    clock_t begin = clock();
    //while(r->t<tmax){
    reb_simulation_integrate(r, tmax);
        //fprintf(f, "%f,%f\n", r->t, ho->y[0]);
    //}
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //printf("\n%e\n", fabs((reb_simulation_energy(r) - e_init) / e_init));
    printf("\n%s %e %e %e\n", grators[integrator], peri, emax, time_spent);

    FILE* f = fopen(title_stats,"a");
    fprintf(f, "%s,%e,%e,%e\n", grators[integrator], peri, emax, time_spent);
    fclose(f);
    reb_simulation_free(r);
}


void heartbeat(struct reb_simulation* r){
  //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
    //if (reb_simulation_output_check(r, 10.*2.*M_PI)){
    //    reb_simulation_output_timing(r, tmax);
    //}

    double e_curr = fabs((reb_simulation_energy(r) - e_init) / e_init);
    if (emax < e_curr){
      emax = e_curr;
    }
/*
    if (reb_simulation_output_check(r, 0.15 * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file

        double err = fabs((reb_simulation_energy(r) - e_init) / e_init);
        if (err > emax){
          emax = err;
        }
        FILE* f = fopen("energy_trace_best.txt","a");

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
    }
*/

}

