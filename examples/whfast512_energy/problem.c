// Usage:  ./rebound <seed> <nsteps> <nsamples>
//   seed=0  -> unperturbed reference; seed>0 -> 1e-15 random perturbation.
#include "rebound.h"
#include "integrator_whfast512.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Solar System initial conditions (NASA Horizons), from examples/whfast512_solar_system.
static double all_ss_pos[9][3] = {
    {-0.008816286905115728, -0.0010954664916791675, 0.0002143249385447027},
    {-0.05942272929227954, -0.46308699693348293, -0.032897989948949075},
    {-0.7276101005375593, 0.006575003332463933, 0.041795901908847084},
    {-0.5789530452882667, -0.8361530119313055, 0.0002611520181901174},
    {-1.45202492400084, 0.827519404194876, 0.052981833432457694},
    {4.492983939852296, 2.0661626247490354, -0.10909246996001629},
    {8.4974210980544, -4.8620394993693585, -0.2537835862373596},
    {12.959111916929283, 14.760785302864473, -0.1130656917933948},
    {29.787987348666505, -2.51460654509393, -0.6347108842010732}
};
static double all_ss_vel[9][3] = {
    {0.00014315746073017681, -0.0004912441820893999, 8.127678560998346e-07},
    {1.2978664284760637, -0.09524541469911743, -0.12677574364801253},
    {-0.019239782390457125, -1.1813975672919448, -0.01509392594251431},
    {0.8098712561282222, -0.5682496529341624, 2.6169897281383047e-05},
    {-0.37436417754222295, -0.6365841544564991, -0.004143932260467942},
    {-0.18818907783656452, 0.41919544951404614, 0.0024710497024424977},
    {0.14292308496870448, 0.2808676923735748, -0.010574288572728728},
    {-0.1734971049470612, 0.14019515029516152, 0.0027683484887051457},
    {0.014142947617173336, 0.18292110872737416, -0.004092845767710294}
};
static double all_ss_mass[9] = {
    0.9999999999950272, 1.6601208254808336e-07, 2.447838287784771e-06,
    3.0404326489511185e-06, 3.2271560828978514e-07, 0.0009547919099366768,
    0.0002858856700231729, 4.366249613200406e-05, 5.151383772628957e-05
};

static double gr_potential(struct reb_simulation* const r){
    const struct reb_particle* const particles = r->particles;
    double C2 = 10065.32 * 10065.32;
    const struct reb_particle source = particles[0];
    const double mu = r->G*source.m;
    const double prefac = 3.*mu*mu/C2;
    double H = 0.;
    for (size_t i=1;i<r->N;i++){
        struct reb_particle pi = particles[i];
        double dx = pi.x - source.x, dy = pi.y - source.y, dz = pi.z - source.z;
        double r2 = dx*dx + dy*dy + dz*dz;
        H -= prefac*pi.m/r2;
    }
    return H;
}

int main(int argc, char** argv){
    setbuf(stdout, NULL);
    long seed   = (argc>1)? atol(argv[1]) : 0;
    long nsteps = (argc>2)? atol(argv[2]) : 200000000L; // total full steps
    int  nsamp  = (argc>3)? atoi(argv[3]) : 100;        // log-spaced samples
    int  gr     = (argc>4)? atoi(argv[4]) : 1;          // GR on/off (debug)
    int  corr   = (argc>5)? atoi(argv[5]) : 17;         // symplectic corrector order
    long conc   = (argc>6)? atol(argv[6]) : 1000;       // concatenate_steps
    const char* integ = (argc>7)? argv[7] : "whfast512"; // integrator (debug)
    double dtdays = (argc>8)? atof(argv[8]) : 6.0;       // timestep in days (debug)

    struct reb_simulation* r = reb_simulation_create();
    r->dt = dtdays/365.25*2*M_PI;       // timestep
    r->G = 1.0;
    r->exact_finish_time = 0;            // integrate in whole steps
    r->force_is_velocity_dependent = 0;
    reb_simulation_set_integrator(r, integ);
    if (strcmp(integ,"whfast512")==0){
        struct reb_integrator_whfast512_state* w = r->integrator.state;
        w->gr_potential = gr;
        w->corrector = corr;
        w->concatenate_steps = conc;
    }

    srand48(seed);
    const double pert = (seed==0) ? 0.0 : 1e-15;
    for (int i=0;i<9;i++){
        struct reb_particle p = {0};
        p.m  = all_ss_mass[i];
        p.x  = all_ss_pos[i][0]*(1.+pert*(2.*drand48()-1.));
        p.y  = all_ss_pos[i][1]*(1.+pert*(2.*drand48()-1.));
        p.z  = all_ss_pos[i][2]*(1.+pert*(2.*drand48()-1.));
        p.vx = all_ss_vel[i][0]*(1.+pert*(2.*drand48()-1.));
        p.vy = all_ss_vel[i][1]*(1.+pert*(2.*drand48()-1.));
        p.vz = all_ss_vel[i][2]*(1.+pert*(2.*drand48()-1.));
        reb_simulation_add(r, p);
    }
    reb_simulation_move_to_com(r);

    const double grE0 = gr ? gr_potential(r) : 0.0;
    const double E0 = reb_simulation_energy(r) + grE0;
    const double dt = r->dt;

    fprintf(stderr, "E0=%.10e  grpot=%.4e  gr=%d seed=%ld\n", E0, grE0, gr, seed);
    printf("# seed=%ld nsteps=%ld dt=%.10e (6 days) gr=%d\n", seed, nsteps, dt, gr);
    printf("# step  t_years  |dE/E|\n");

    const double s0 = 1000.0;
    long last = 0;
    for (int s=1;s<=nsamp;s++){
        double frac = (double)s/(double)nsamp;
        long target = (long)(s0 * pow((double)nsteps/s0, frac));
        if (target<=last) continue;           // keep strictly increasing
        last = target;
        reb_simulation_integrate(r, target*dt);
        reb_simulation_synchronize(r);
        double E = reb_simulation_energy(r) + (gr ? gr_potential(r) : 0.0);
        double tyears = r->t/(2.*M_PI);
        printf("%ld %.6e %.8e\n", target, tyears, fabs((E-E0)/E0));
    }
    reb_simulation_free(r);
    return 0;
}
