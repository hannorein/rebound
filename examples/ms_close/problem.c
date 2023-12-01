/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
//#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "integrator_trace.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double err;
double emax = 0.0;

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_create_simulation();
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
    //r->ri_trace.S_peri = reb_integrator_trace_switch_fdot_peri;
    //r->ri_trace.S_peri = reb_integrator_trace_peri_switch_default;
    //r->ri_trace.vfac_p = 32.0;

    r->dt = 0.15*2.*M_PI;
    r->heartbeat = heartbeat;


    reb_add(r, star);
    reb_add_fmt(r, "m a e", jm, ja, je);

    //r->dt = orb.P / 15.12345;
    reb_add_fmt(r, "primary m a e inc pomega f", star, sm, sa, se, si, M_PI/2, 0.);
    struct reb_orbit orb = reb_tools_particle_to_orbit(r->G, r->particles[2], r->particles[0]);

    reb_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    r->heartbeat  = heartbeat;
    e_init = reb_tools_energy(r);
    system("rm -rf energy_fdot.txt");
    //FILE* f = fopen("energy_wh.txt","w");
    //fprintf(f, "t,sf,stauf,jf,jtauf\n");

    reb_integrate(r, 300*2*M_PI*29.4);
    //reb_integrate(r, orb.P * 3.);
    //reb_integrate(r,2. *M_PI*11.86);
    //reb_integrate(r, 1277.);
    //err = reb_tools_energy(r);
    printf("%e\n", emax);
    reb_free_simulation(r);
}


void heartbeat(struct reb_simulation* r){
    double err = fabs((reb_tools_energy(r) - e_init) / e_init);
    if (err > emax){
      emax = err;
    }
  //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
    //if (reb_output_check(r, 10.*2.*M_PI)){
    //    reb_output_timing(r, 0);
    //}

    //if (reb_output_check(r, 10. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file
        FILE* f = fopen("energy_fdot.txt","a");

        // rotate whole simulation to rotating frame
        //reb_simulation_irotate(r, r1);

        struct reb_particle* sun = &r->particles[0];
        struct reb_particle* jup = &r->particles[1];
        struct reb_particle* sat = &r->particles[2];

        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, *sat, *sun);
        struct reb_orbit oj = reb_tools_particle_to_orbit(r->G, *jup, *sun);

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
        double e = (reb_tools_energy(r) - e_init) / e_init;
        fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %f %f %f %f\n",r->t, e, r->particles[0].x, r->particles[0].y, r->particles[0].z, r->particles[1].x, r->particles[1].y, r->particles[1].z, r->particles[2].x, r->particles[2].y, r->particles[2].z, o.e, o.inc, dvx/vcirc, dvxj/vcircj);

        //double e = reb_tools_energy(r);
        //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
        fclose(f);

        //reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);
    //}

}
