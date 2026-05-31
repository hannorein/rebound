#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "rebound.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline double rand_unit(uint64_t* s0, uint64_t* s1) {
    uint64_t result = *s0 + *s1;
    uint64_t s1x = *s1 ^ *s0;
    *s0 = rotl(*s0, 24) ^ s1x ^ (s1x << 16);
    *s1 = rotl(s1x, 37);
    return (result >> 11) * (1.0 / (double)(1ULL << 53));
}

// Vanishing planet mass: asm512 rejects test particles (N_active != N).
static const double TEST_MASS = 1e-30;

// One cell per simulation (N=2). asm512's mat8_mul3 couples all 8 SIMD
// lanes, so cells batched together can contaminate each other when one
// lane diverges. Using N=2 keeps only lane 0 active (mask=0x01) and makes
// every cell evaluate independently.
static double run_cell_asm512(const char* integrator_name,
                              double a0, double e0,
                              double dt, double tmax, uint64_t seed) {
    uint64_t s0 = seed ^ 0x9E3779B97F4A7C15ULL;
    uint64_t s1 = (seed * 0x6A5D39EAE12657AAULL) ^ 0xBF58476D1CE4E5B9ULL;
    (void)rand_unit(&s0, &s1);
    double f0 = 2.0 * M_PI * rand_unit(&s0, &s1);

    struct reb_simulation* r = reb_simulation_create();
    r->G = 1.0;
    r->dt = dt;
    r->exact_finish_time = 0;
    reb_simulation_set_integrator(r, integrator_name);
    struct reb_integrator_asm512_state* st = r->integrator.state;
    st->gr_potential = 0;
    st->N_systems    = 1;

    struct reb_particle star = {0};
    star.m = 1.0;
    reb_simulation_add(r, star);

    struct reb_particle tp = reb_particle_from_orbit(
        r->G, star, TEST_MASS, a0, e0, 0.0, 0.0, 0.0, f0);
    tp.m = TEST_MASS;
    reb_simulation_add(r, tp);

    reb_simulation_integrate(r, tmax);
    reb_simulation_synchronize(r);

    double da_over_a = NAN;
    if (r->N == 2) {
        struct reb_orbit o = reb_orbit_from_particle(
            r->G, r->particles[1], r->particles[0]);
        if (isfinite(o.a) && o.a > 0.0) {
            da_over_a = (o.a - a0) / a0;
        }
    }
    reb_simulation_free(r);
    return da_over_a;
}

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [--integrator NAME] [--na N] [--ne N]\n"
        "          [--amin x] [--amax x] [--emin x] [--emax x]\n"
        "          [--dt-days x] [--t-years x] [--seed N] [--output FILE]\n"
        "Integrators: asm512 | asm512_opt\n", prog);
}

int main(int argc, char* argv[]) {
    int na = 200, ne = 200;
    double amin = 0.05, amax = 1.0;
    double emin = 0.0,  emax = 0.95;
    double dt_days = 5.0;
    double t_years = 10.0;
    uint64_t seed = 42;
    const char* output = "kepler_stability.csv";
    const char* integ_name = "asm512";

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--na") && i + 1 < argc) na = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--ne") && i + 1 < argc) ne = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--amin") && i + 1 < argc) amin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--amax") && i + 1 < argc) amax = atof(argv[++i]);
        else if (!strcmp(argv[i], "--emin") && i + 1 < argc) emin = atof(argv[++i]);
        else if (!strcmp(argv[i], "--emax") && i + 1 < argc) emax = atof(argv[++i]);
        else if (!strcmp(argv[i], "--dt-days") && i + 1 < argc) dt_days = atof(argv[++i]);
        else if (!strcmp(argv[i], "--t-years") && i + 1 < argc) t_years = atof(argv[++i]);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = (uint64_t)atoll(argv[++i]);
        else if (!strcmp(argv[i], "--output") && i + 1 < argc) output = argv[++i];
        else if (!strcmp(argv[i], "--integrator") && i + 1 < argc) {
            integ_name = argv[++i];
            if (strcmp(integ_name, "asm512") && strcmp(integ_name, "asm512_opt")) {
                fprintf(stderr, "Unknown integrator: %s\n", integ_name); return 1;
            }
        }
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
        else { fprintf(stderr, "Unknown arg: %s\n", argv[i]); usage(argv[0]); return 1; }
    }

    if (na < 2 || ne < 2) { fprintf(stderr, "na and ne must be >= 2\n"); return 1; }

    const double dt   = dt_days / 365.25 * 2.0 * M_PI;
    const double tmax = t_years           * 2.0 * M_PI;

    fprintf(stderr,
        "Sweep: integrator=%s  na=%d ne=%d  a in [%.4f, %.4f]  e in [%.4f, %.4f]\n"
        "       dt = %.3f days (= %.6e code)  tmax = %.2f yr (= %.6e code)\n"
        "       seed = %llu  output = %s\n",
        integ_name, na, ne, amin, amax, emin, emax,
        dt_days, dt, t_years, tmax,
        (unsigned long long)seed, output);
#ifdef _OPENMP
    fprintf(stderr, "       OpenMP threads = %d\n", omp_get_max_threads());
#endif

    const int Ntot = na * ne;
    double* da = (double*)malloc((size_t)Ntot * sizeof(double));
    if (!da) { fprintf(stderr, "malloc failed\n"); return 1; }

    const double t_start = (double)clock() / CLOCKS_PER_SEC;
    double t_wall0 = 0.0;
#ifdef _OPENMP
    t_wall0 = omp_get_wtime();
#endif

    int done = 0;
    #pragma omp parallel for schedule(dynamic, 8)
    for (int idx = 0; idx < Ntot; idx++) {
        int ia = idx / ne;
        int ie = idx % ne;
        double a0 = amin + (amax - amin) * (double)ia / (double)(na - 1);
        double e0 = emin + (emax - emin) * (double)ie / (double)(ne - 1);
        uint64_t cell_seed = seed * 1000003ULL + (uint64_t)idx;
        da[idx] = run_cell_asm512(integ_name, a0, e0, dt, tmax, cell_seed);

        #pragma omp atomic
        done++;
        if ((done & 1023) == 0) {
            fprintf(stderr, "\r  progress: %6d / %d  (%.1f%%)",
                    done, Ntot, 100.0 * done / Ntot);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\r  progress: %6d / %d  (100.0%%)\n", Ntot, Ntot);

    double t_wall1 = 0.0;
#ifdef _OPENMP
    t_wall1 = omp_get_wtime();
#endif
    fprintf(stderr, "  wall = %.2f s  (cpu %.2f s)\n",
            t_wall1 - t_wall0, (double)clock() / CLOCKS_PER_SEC - t_start);

    FILE* fp = fopen(output, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", output); free(da); return 1; }
    fprintf(fp, "# na=%d ne=%d amin=%.10g amax=%.10g emin=%.10g emax=%.10g\n",
            na, ne, amin, amax, emin, emax);
    fprintf(fp, "# dt_days=%.10g t_years=%.10g seed=%llu integrator=%s\n",
            dt_days, t_years, (unsigned long long)seed, integ_name);
    fprintf(fp, "a_init,e_init,da_over_a\n");
    for (int ia = 0; ia < na; ia++) {
        double a0 = amin + (amax - amin) * (double)ia / (double)(na - 1);
        for (int ie = 0; ie < ne; ie++) {
            double e0 = emin + (emax - emin) * (double)ie / (double)(ne - 1);
            double v  = da[ia * ne + ie];
            if (isnan(v)) {
                fprintf(fp, "%.10g,%.10g,nan\n", a0, e0);
            } else {
                fprintf(fp, "%.10g,%.10g,%.16e\n", a0, e0, v);
            }
        }
    }
    fclose(fp);

    int nfin = 0, nnan = 0;
    for (int i = 0; i < Ntot; i++) { if (isnan(da[i])) nnan++; else nfin++; }
    fprintf(stderr, "  wrote %s : %d finite, %d NaN\n", output, nfin, nnan);

    free(da);
    return 0;
}
