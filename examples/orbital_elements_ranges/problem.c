/**
 * Check the definition ranges of the angular orbital elements.
 */
#include "rebound.h"
#include <assert.h>
#include <stdio.h>
#include <math.h>

struct angle_range {
    double min;
    double max;
};

enum angle_index {
    INC,
    OMEGA,
    ARG_PERI,
    POMEGA,
    THETA,
    TRUE_ANOM,
    MEAN_ANOM,
    ECC_ANOM,
    MEAN_LONG,
    ANGLE_COUNT
};

static void update_range(struct angle_range* range, double value){
    if (value < range->min) range->min = value;
    if (value > range->max) range->max = value;
}

static void update_orbit_ranges(struct angle_range* ranges, struct reb_orbit orbit){
    update_range(&ranges[INC], orbit.inc);
    update_range(&ranges[OMEGA], orbit.Omega);
    update_range(&ranges[ARG_PERI], orbit.omega);
    update_range(&ranges[POMEGA], orbit.pomega);
    update_range(&ranges[THETA], orbit.theta);
    update_range(&ranges[TRUE_ANOM], orbit.f);
    update_range(&ranges[MEAN_ANOM], orbit.M);
    update_range(&ranges[ECC_ANOM], reb_M_to_E(orbit.e, orbit.M));
    update_range(&ranges[MEAN_LONG], orbit.l);
}

int main(void){
    // Init range for each angle
    struct angle_range ranges[ANGLE_COUNT];
    for (unsigned int i = 0; i < ANGLE_COUNT; i++){
        ranges[i].min = INFINITY;
        ranges[i].max = -INFINITY;
    }

    // Init simulation with primary
    struct reb_particle primary = {0};
    primary.m = 1.;
    struct reb_simulation* simulation = reb_simulation_create();
    simulation->rand_seed = 0;

    // Define hand-picked test cases
    const double pi = M_PI;
    const double pi2 = 2.*M_PI;
    const double tol = 1.e-6;
    const double hand_picked[][4] = {
        {0., 0., 0., 0.},
        {pi, pi, 0., 0.},
        {pi/4., pi, pi, 0.},
        {pi/4., pi - tol, pi, 0.},
        {3.*pi/4., pi + tol, pi, 0.},
        {pi/2., pi + tol, 0., 0.},
        {pi/2., -pi2 + tol, -pi2 + tol, -pi2 + tol},
        {pi/2., pi2 - tol, pi2 - tol, pi2 - tol}
    };

    // Update ranges for hand-picked test cases
    for (unsigned int i = 0; i < sizeof(hand_picked)/sizeof(hand_picked[0]); i++){
        struct reb_particle particle = reb_particle_from_orbit(
            1., primary, 0., 1., 0.5,
            hand_picked[i][0], hand_picked[i][1],
            hand_picked[i][2], hand_picked[i][3]);
        update_orbit_ranges(ranges, reb_orbit_from_particle(1., particle, primary));
    }

    // Update ranges for random test cases
    const unsigned int N = 10000;
    const double angle_min = -10.*M_PI;
    const double angle_max = 10.*M_PI;
    for (unsigned int i = 0; i < N; i++){
        const double e = reb_random_uniform(simulation, 0., 0.9);
        const double inc = reb_random_uniform(simulation, angle_min, angle_max);
        const double Omega = reb_random_uniform(simulation, angle_min, angle_max);
        const double omega = reb_random_uniform(simulation, angle_min, angle_max);
        const double f = reb_random_uniform(simulation, angle_min, angle_max);
        struct reb_particle particle = reb_particle_from_orbit(
            1., primary, 0., 1., e, inc, Omega, omega, f);
        update_orbit_ranges(ranges, reb_orbit_from_particle(1., particle, primary));
    }

    // Print results
    printf("Orbital elements ranges after test cases:\n");
    printf("inc      min=%.17g  max=%.17g\n", ranges[INC].min, ranges[INC].max);
    printf("Omega    min=%.17g  max=%.17g\n", ranges[OMEGA].min, ranges[OMEGA].max);
    printf("omega    min=%.17g  max=%.17g\n", ranges[ARG_PERI].min, ranges[ARG_PERI].max);
    printf("pomega   min=%.17g  max=%.17g\n", ranges[POMEGA].min, ranges[POMEGA].max);
    printf("theta    min=%.17g  max=%.17g\n", ranges[THETA].min, ranges[THETA].max);
    printf("f        min=%.17g  max=%.17g\n", ranges[TRUE_ANOM].min, ranges[TRUE_ANOM].max);
    printf("M        min=%.17g  max=%.17g\n", ranges[MEAN_ANOM].min, ranges[MEAN_ANOM].max);
    printf("E        min=%.17g  max=%.17g\n", ranges[ECC_ANOM].min, ranges[ECC_ANOM].max);
    printf("l        min=%.17g  max=%.17g\n", ranges[MEAN_LONG].min, ranges[MEAN_LONG].max);

    // Assert expected ranges
    printf("\nAsserting expected ranges for orbital elements...\n");
    assert(ranges[INC].min >= 0. && ranges[INC].max <= pi);
    assert(ranges[OMEGA].min > -pi && ranges[OMEGA].max <= pi);
    assert(ranges[ARG_PERI].min >= 0. && ranges[ARG_PERI].max < pi2);
    assert(ranges[POMEGA].min > -pi2 && ranges[POMEGA].max <= pi2);
    assert(ranges[THETA].min >= 0. && ranges[THETA].max < pi2);
    assert(ranges[TRUE_ANOM].min >= 0. && ranges[TRUE_ANOM].max < pi2);
    assert(ranges[MEAN_ANOM].min >= 0. && ranges[MEAN_ANOM].max < pi2);
    assert(ranges[ECC_ANOM].min >= 0. && ranges[ECC_ANOM].max < pi2);
    assert(ranges[MEAN_LONG].min >= 0. && ranges[MEAN_LONG].max < pi2);

    // End of test
    printf("All orbital elements are within the expected ranges.\n");
    reb_simulation_free(simulation);
}