#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);

static const char* PLANET[8] = {
    "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune",
};

int main(int argc, char* argv[]) {
    double t_end_yr = (argc > 1) ? atof(argv[1]) : 1000.0;

    double dt_days = (argc > 2) ? atof(argv[2]) : 5.0;
    struct reb_simulation* r = reb_simulation_create();
    r->dt = dt_days / 365.25 * 2 * M_PI;
    r->G  = 1.0;
    r->exact_finish_time = 0;
    reb_simulation_add_fmt(r, "solarsystem");

    reb_simulation_set_integrator(r, "asm512");
    struct reb_integrator_asm512_state* st = r->integrator.state;
    st->gr_potential = 1;
    st->concatenate_steps = 1000000;
    st->corrector = 17;

    long n_steps = (long)(t_end_yr * 365.25 / dt_days);
    reb_simulation_integrate(r, t_end_yr * 2 * M_PI);

    printf("# t_end = %.3e yr   dt = %.2f d   total_steps ~ %ld\n",
           t_end_yr, dt_days, n_steps);
    printf("# lane  planet      newton+bisection_count   per_step\n");
    uint64_t total = 0;
    uint64_t lane_counts[8];
    for (int k = 0; k < 8; k++) {
        lane_counts[k] = reb_asm512_counter(r, k);
        total += lane_counts[k];
    }
    for (int k = 0; k < 8; k++) {
        printf("%6d  %-8s    %18" PRIu64 "   %10.4f\n",
               k, PLANET[k], lane_counts[k],
               (double)lane_counts[k] / (double)n_steps);
    }
    printf("# total iter-particle-events: %" PRIu64 "\n", total);
    printf("# average per step (sum over 8 lanes): %.4f\n",
           (double)total / (double)n_steps);

    reb_simulation_free(r);
    return 0;
}
