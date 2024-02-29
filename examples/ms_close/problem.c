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
double etot;
double tmax = 300. * 2. * M_PI * 29.4;

int ind;
char title[100] = "test.txt";
char title_remove[100] = "rm -rf test.txt";

int peri_switch = 3;
int ts;
int eta;
double tracker;

int main(int argc, char* argv[]){

    double eccs[15] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999};
    double etas[10] = {0.1, 0.5, 1.0, 2.0, 3.0};
    double timesteps[5] = {0.15 * 2 * M_PI, 11.8 / 30., 11.8 / 20.};

    if (argc == 3){
      ind = atoi(argv[2]);
      ts = atoi(argv[1]);
    }

    if (argc == 4){
      eta = atoi(argv[3]);
      ind = atoi(argv[2]);
      ts = atoi(argv[1]);
    }

    // Initialize masses
    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle star = {0};
    star.m = 1;

    // Jupiter
    double jm = 9.55e-4;
    double ja = 5.2;
    double je = 0.0;
    double jmse = 1.-je;

    // Saturn
    double sm = 2.857e-4;
    double sa = 9.58;
    double se = eccs[ind];
    double omse = 1. - se;
    double si = M_PI / 2.;

    reb_simulation_add(r, star);
    reb_simulation_add_fmt(r, "m a e", jm, ja, je);

    reb_simulation_add_fmt(r, "primary m a e inc omega f", star, sm, sa, se, si, M_PI/2.,M_PI);

    r->integrator = REB_INTEGRATOR_TRACE;
    r->dt = timesteps[ts];
    r->ri_trace.peri_mode = REB_TRACE_PERI_PARTIAL_BS;

    switch (peri_switch){
      case 1:
        r->ri_trace.S_peri = reb_integrator_trace_switch_peri_pham2024;
        r->ri_trace.peri_crit_eta = etas[eta];
        tracker = r->ri_trace.peri_crit_eta;
        //printf("case 1\n");
        break;
      case 2:
        r->ri_trace.S_peri = reb_integrator_trace_switch_peri_default;
        tracker = r->ri_trace.peri_crit_fdot;
        //printf("case 2\n");
        break;
      case 3:
        r->ri_trace.S_peri = reb_integrator_trace_switch_peri_distance;
        r->ri_trace.peri_crit_distance = 2.;
        tracker = r->ri_trace.peri_crit_distance;
        break;
    }

    r->heartbeat = heartbeat;

    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_simulation_energy(r);
    etot = 0.0;
    system(title_remove);

    clock_t begin = clock();
    reb_simulation_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    int steps_done = r->steps_done;

    FILE* f = fopen("switch_stats.txt","a");
    fprintf(f, "%e,%f,%f,%f,%f\n", etot/steps_done, se, ts,r->dt,tracker);
    fclose(f);
    reb_simulation_free(r);
}


void heartbeat(struct reb_simulation* r){
  //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
    //if (reb_simulation_output_check(r, 10.*2.*M_PI)){
    //    reb_simulation_output_timing(r, tmax);
    //}
    //printf("etot: %e %e\n", etot, fabs((reb_simulation_energy(r) - e_init) / e_init));

    etot += fabs((reb_simulation_energy(r) - e_init) / e_init);
/*
    //if (reb_simulation_output_check(r, 10. * 2.*M_PI)){
        // Once per 4 days, output the relative energy error to a text file

        //FILE* f = fopen("energy_trace_partial_bs.txt","a");

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

        //const double vcirc = sqrt(r->G * r->particles[0].m / o.a);
        //const double vcircj = sqrt(r->G * r->particles[0].m / oj.a);
        double e = (reb_simulation_energy(r) - e_init) / e_init;

        //double e = reb_tools_energy(r);
        //fprintf(f,"%e %e\n",r->t, (e - e_init) / e_init);
        FILE* f = fopen(title,"a");
        fprintf(f,"%f,%0.20e,%f,%f\n",r->t, e, o.e, o.inc);
        fclose(f);

        //reb_integrator_synchronize(r);
        //printf("\n Jacobi: %e\n", e);
    //}

/*
    struct reb_orbit orb = reb_orbit_from_particle(r->G, r->particles[1], r->particles[0]);

    struct reb_particle s = r->particles[0];
    struct reb_particle p = r->particles[2];

    int j = 2;
    double GM = r->G*r->particles[0].m; // Not sure if this is the right mass to use.
    double x = r->particles[j].x - r->particles[0].x;
    double y = r->particles[j].y - r->particles[0].y;
    double z = r->particles[j].z - r->particles[0].z;
    double d2 = x*x + y*y + z*z;

    // first derivative
    double dx = r->particles[j].vx - r->particles[0].vx;
    double dy = r->particles[j].vy - r->particles[0].vy;
    double dz = r->particles[j].vz - r->particles[0].vz;

    // second derivative
    double prefact2 = -GM/sqrt(d2*d2*d2);
    double ddx = prefact2*x;
    double ddy = prefact2*y;
    double ddz = prefact2*z;
    double dd = sqrt(ddx*ddx + ddy*ddy + ddz*ddz);

    // third derivative
    double prefact3 = GM/sqrt(d2*d2*d2*d2*d2);
    double dddx = prefact3*(-dx*(y*y+z*z) + 2.*x*x*dx+3.*x*(y*dy+z*dz));
    double dddy = prefact3*(-dy*(x*x+z*z) + 2.*y*y*dy+3.*y*(x*dx+z*dz));
    double dddz = prefact3*(-dz*(x*x+y*y) + 2.*z*z*dz+3.*z*(x*dx+y*dy));
    double ddd = sqrt(dddx*dddx + dddy*dddy + dddz*dddz);

    // fourth derivative
    double prefact4 = GM/sqrt(d2*d2*d2*d2*d2*d2*d2);
    double ddddx = prefact4* (d2 * (-ddx*(y*y+z*z) + 2.*x*x*ddx + dx*(y*dy + z*dz) + x*(4.*dx*dx + 3.*(y*ddy + dy*dy + z*ddz + dz*dz ))) - 5.*(x*dx+y*dy+z*dz)*(-dx*(y*y+z*z)+2.*x*x*dx + 3.*x*(y*dy+z*dz)));
    double ddddy = prefact4* (d2 * (-ddy*(x*x+z*z) + 2.*y*y*ddy + dy*(x*dx + z*dz) + y*(4.*dy*dy + 3.*(x*ddx + dx*dx + z*ddz + dz*dz ))) - 5.*(y*dy+x*dx+z*dz)*(-dy*(x*x+z*z)+2.*y*y*dy + 3.*y*(x*dx+z*dz)));
    double ddddz = prefact4* (d2 * (-ddz*(y*y+x*x) + 2.*z*z*ddz + dz*(y*dy + x*dx) + z*(4.*dz*dz + 3.*(y*ddy + dy*dy + x*ddx + dx*dx ))) - 5.*(z*dz+y*dy+x*dx)*(-dz*(y*y+x*x)+2.*z*z*dz + 3.*z*(y*dy+x*dx)));
    double dddd = sqrt(ddddx*ddddx + ddddy*ddddy + ddddz*ddddz);

    double tau_prs23 = sqrt(2.*dd*dd/(ddd*ddd+dd*dddd)); // Eq 16
    //double eta = 0.1; // Requires experimentation

    const double hx = (y*dz - z*dy);  // specific angular momentum vector
    const double hy = (z*dx - x*dz);
    const double hz = (x*dy - y*dx);
    const double h2 = hx*hx + hy*hy + hz*hz;

    // This only works for bound orbits!
    const double fdot2 = h2 / (d2*d2);
    const double peff2 = (4 * M_PI * M_PI) / fdot2; // effective period squared

    FILE* f = fopen(title,"a");
    fprintf(f,"%f,%e,%f,%f\n",r->t, tau_prs23,sqrt(d2),sqrt(peff2));
    fclose(f);
*/

}
