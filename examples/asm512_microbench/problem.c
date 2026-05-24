#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

static struct reb_simulation* setup_solar(void){
    struct reb_simulation* r = reb_simulation_create();
    r->dt = 5.0/365.25*2*M_PI;
    r->G  = 1.;
    r->exact_finish_time = 0;
    reb_simulation_add_fmt(r, "solarsystem");
    reb_simulation_set_integrator(r, "asm512");
    struct reb_integrator_asm512_state* asm512 = r->integrator.state;
    asm512->gr_potential      = 0;
    asm512->concatenate_steps = 1000000;
    asm512->corrector         = 17;
    return r;
}

static double wall_seconds(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec*1e-9;
}

static int cmp_double(const void* a, const void* b){
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

int main(int argc, char* argv[]){
    setvbuf(stdout, NULL, _IONBF, 0);
    const char* tag    = (argc > 1) ? argv[1] : "asm512";
    double tmax_yr     = (argc > 2) ? atof(argv[2]) : 5e3;
    int trials         = (argc > 3) ? atoi(argv[3]) : 11;
    const char* csv    = (argc > 4) ? argv[4] : "results.csv";
    if (trials < 1)  trials = 1;
    if (trials > 31) trials = 31;

    double times[32];
    double e_err = 0;

    {
        struct reb_simulation* r = setup_solar();
        reb_simulation_integrate(r, tmax_yr*2*M_PI);
        reb_simulation_free(r);
    }

    for (int t = 0; t < trials; t++){
        struct reb_simulation* r = setup_solar();
        double E0 = reb_simulation_energy(r);
        double t0 = wall_seconds();
        reb_simulation_integrate(r, tmax_yr*2*M_PI);
        double t1 = wall_seconds();
        times[t] = t1 - t0;
        if (t == trials-1){
            double E1 = reb_simulation_energy(r);
            e_err = fabs((E0 - E1)/E0);
        }
        reb_simulation_free(r);
    }

    qsort(times, trials, sizeof(double), cmp_double);
    double t_min = times[0];
    double t_med = times[trials/2];
    double t_max = times[trials-1];

    printf("%-24s  trials=%2d  min=%.4fs  med=%.4fs  max=%.4fs  E_rel=%.3e\n",
           tag, trials, t_min, t_med, t_max, e_err);

    // Append to CSV (header if new)
    int need_header = 0;
    FILE* f = fopen(csv, "r");
    if (!f) { need_header = 1; }
    else fclose(f);
    f = fopen(csv, "a");
    if (!f){ perror(csv); return 1; }
    if (need_header) fprintf(f, "tag,trials,t_min,t_med,t_max,e_rel_err,tmax_yr\n");
    fprintf(f, "%s,%d,%.6f,%.6f,%.6f,%.6e,%.6e\n",
            tag, trials, t_min, t_med, t_max, e_err, tmax_yr);
    fclose(f);

    return 0;
}
