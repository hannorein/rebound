#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);

#define NCFG 4
static const char* configs[NCFG]      = {"asm512", "asm512_mom", "asm512_fused", "asm512_opt"};
static const int   instrumented[NCFG] = {0, 1, 1, 1};

static struct reb_simulation* setup_sim_solar(const char* integ_name){
    struct reb_simulation* r = reb_simulation_create();
    r->dt = 5.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
    reb_simulation_add_fmt(r, "solarsystem");
    reb_simulation_set_integrator(r, integ_name);
    struct reb_integrator_asm512_state* asm512 = r->integrator.state;
    asm512->gr_potential    = 0;
    asm512->concatenate_steps = 1e6;
    asm512->corrector       = 17;
    return r;
}

// 8 planets, small a, low e, dt sized so dt/T is in the regime where Newton
// must do 1-2 iters (no bisection). e=0.05 keeps orbits smooth so momentum's
// linear extrapolation should be a much better guess than dt/r.
static struct reb_simulation* setup_sim_stress(const char* integ_name){
    struct reb_simulation* r = reb_simulation_create();
    r->dt = 25.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
    r->rand_seed = 42;
    reb_simulation_add_fmt(r, "m", 1.0);
    const double a_vals[8] = {0.30, 0.33, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56};
    for (int i = 0; i < 8; i++){
        reb_simulation_add_fmt(r, "a e uniform(f)", a_vals[i], 0.05);
    }
    reb_simulation_set_integrator(r, integ_name);
    struct reb_integrator_asm512_state* asm512 = r->integrator.state;
    asm512->gr_potential    = 0;
    asm512->concatenate_steps = 1e6;
    asm512->corrector       = 17;
    return r;
}

typedef struct reb_simulation* (*setup_fn)(const char*);
static setup_fn pick_setup(const char* name){
    if (strcmp(name, "stress") == 0) return setup_sim_stress;
    return setup_sim_solar;
}

static double wall_seconds(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec*1e-9;
}

static uint64_t total_counter(struct reb_simulation* r, int n_active_lanes){
    uint64_t s = 0;
    for (int i=0; i<n_active_lanes; i++) s += reb_asm512_counter(r, i);
    return s;
}

static int cmp_double(const void* a, const void* b){
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

static void run_bench(const char* setup_name, double tmax_yr, int trials){
    const int n_active = 8;
    setup_fn setup = pick_setup(setup_name);
    printf("# setup=%s, corrector=17, no GR, tmax=%.2g yr, trials=%d\n", setup_name, tmax_yr, trials);
    printf("# counter = extra Newton iters (above first convergence) summed over all steps + lanes\n");
    printf("%-15s %12s %14s %18s %14s\n", "config", "walltime_s", "E_rel_err", "counter_total", "ctr/step/lane");
    for (int c = 0; c < NCFG; c++){
        double times[16];
        double e_err = 0;
        uint64_t ctot = 0;
        double n_steps_per_lane = 0;
        for (int t = 0; t < trials; t++){
            struct reb_simulation* r = setup(configs[c]);
            double E0 = reb_simulation_energy(r);
            double t0 = wall_seconds();
            reb_simulation_integrate(r, tmax_yr*2*M_PI);
            double t1 = wall_seconds();
            times[t] = t1 - t0;
            if (t == trials-1){
                double E1 = reb_simulation_energy(r);
                e_err = fabs((E0 - E1)/E0);
                ctot = instrumented[c] ? total_counter(r, n_active) : 0;
                n_steps_per_lane = r->t / r->dt;
            }
            reb_simulation_free(r);
        }
        qsort(times, trials, sizeof(double), cmp_double);
        double median = times[trials/2];
        double ctr_per = instrumented[c] ? ((double)ctot)/(n_active*n_steps_per_lane) : NAN;
        printf("%-15s %12.4f %14.4e %18lu %14.6f\n",
               configs[c], median, e_err, (unsigned long)ctot, ctr_per);
    }
}

static void run_chunks(const char* setup_name, double tmax_yr, int n_chunks, const char* csv_path){
    const int n_active = 8;
    setup_fn setup = pick_setup(setup_name);
    FILE* f = fopen(csv_path, "w");
    if (!f){ perror(csv_path); exit(1); }
    fprintf(f, "config,t_yr,wall_s_chunk,counter_chunk,ctr_per_step_per_lane\n");
    const double dt_chunk = tmax_yr * 2 * M_PI / n_chunks;
    for (int c = 0; c < NCFG; c++){
        struct reb_simulation* r = setup(configs[c]);
        struct reb_integrator_asm512_state* asm512 = r->integrator.state;
        unsigned int steps_per_chunk = (unsigned int)(dt_chunk / r->dt);
        asm512->concatenate_steps = (steps_per_chunk > 16) ? (steps_per_chunk / 16) : 1;
        uint64_t prev_ctr = 0;
        double prev_t = 0;
        for (int k = 1; k <= n_chunks; k++){
            double t_target = k * dt_chunk;
            double t0 = wall_seconds();
            reb_simulation_integrate(r, t_target);
            double t1 = wall_seconds();
            uint64_t ctr = instrumented[c] ? total_counter(r, n_active) : 0;
            uint64_t dctr = ctr - prev_ctr;
            prev_ctr = ctr;
            double actual_dt = r->t - prev_t;
            prev_t = r->t;
            double steps_in_chunk = actual_dt / r->dt;
            double ctr_per = (instrumented[c] && steps_in_chunk > 0) ? ((double)dctr)/(n_active*steps_in_chunk) : NAN;
            fprintf(f, "%s,%.6e,%.6f,%lu,%.6f\n",
                    configs[c], r->t/(2*M_PI), t1-t0, (unsigned long)dctr, ctr_per);
        }
        reb_simulation_free(r);
        fprintf(stderr, "  %s done\n", configs[c]);
    }
    fclose(f);
    fprintf(stderr, "wrote %s\n", csv_path);
}

int main(int argc, char* argv[]){
    setvbuf(stdout, NULL, _IONBF, 0);
    const char* mode  = (argc > 1) ? argv[1] : "bench";
    const char* setup = (argc > 2) ? argv[2] : "solar";
    double tmax_yr    = (argc > 3) ? atof(argv[3]) : 1e4;
    if (strcmp(mode, "bench") == 0){
        int trials = (argc > 4) ? atoi(argv[4]) : 5;
        if (trials < 1 || trials > 15) trials = 5;
        run_bench(setup, tmax_yr, trials);
    } else if (strcmp(mode, "chunks") == 0){
        int n_chunks    = (argc > 4) ? atoi(argv[4]) : 40;
        const char* csv = (argc > 5) ? argv[5] : "chunks.csv";
        run_chunks(setup, tmax_yr, n_chunks, csv);
    } else {
        fprintf(stderr, "usage: %s bench  [solar|stress] [tmax_yr] [trials]\n", argv[0]);
        fprintf(stderr, "       %s chunks [solar|stress] [tmax_yr] [n_chunks] [csv]\n", argv[0]);
        return 1;
    }
    return 0;
}
