#include "rebound.h"
#include "integrator_whfast.h"
#include "transformations.h"
#include "integrator_symba.h"
#include "assert.h"
#include <math.h>



static void* create_symba_state(void) {
    struct reb_integrator_symba_state *state = malloc(sizeof(struct reb_integrator_symba_state));
    if(state != NULL) {
        state->levels = 3;
        state->r_outer = 0.15;
        state->ratio = 2.08;
    }
    return state;

}

static void fail_step(struct reb_simulation* r, const char* message){
    r->status = REB_STATUS_GENERIC_ERROR;
    reb_simulation_error(r, message);
}

static void sun_drift(const struct reb_simulation* r, struct reb_particle* dh_particles, double dt){
    struct reb_integrator_whfast_state jump = {.coordinates = REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC, .p_jh = dh_particles};
    reb_integrator_whfast_jump_step(r, &jump, dt);
}

static double weight_changing_function(double distance, double r_outer, double r_inner) {
    assert(r_outer > r_inner && r_inner > 0.0);
    assert(distance >= 0.0);
    if (distance <= r_inner) {
        return 0.0;
    }
    if (distance >= r_outer) {
        return 1.0;
    }

    const double w = (r_outer - distance) / (r_outer - r_inner);
    return 2 * w * w * w - 3 * w * w + 1;

}

static double level_weight(double distance, const double *radii, unsigned int level,
        unsigned int max_level){
    assert(level <= max_level);
    if (max_level == 0) {
        return 1.0;
    }
    if(level == 0){
        return weight_changing_function(distance, radii[0], radii[1]);
    }
    if(level == max_level){
        return 1 - weight_changing_function(distance, radii[max_level - 1], radii[max_level]);
    }
    else{
        return weight_changing_function(distance, radii[level], radii[level + 1]) - weight_changing_function(distance, radii[level - 1], radii[level]);
    }
}



struct SyMBA_Pair{
    size_t i;
    size_t j;
};



static void level_wise_interaction_kick(const struct reb_simulation* r, struct reb_particle* dh_particles, const struct SyMBA_Pair* pairs, size_t num_pairs,const double* radii,unsigned int level, unsigned int max_level,double dt){
    for(size_t k = 0; k < num_pairs; k++) {
        struct reb_particle* a = &dh_particles[pairs[k].i];
        struct reb_particle* b = &dh_particles[pairs[k].j];
        const double dx = b->x - a->x;
        const double dy = b->y - a->y;
        const double dz = b->z - a->z;
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double distance = sqrt(r2);
        const double w = level_weight(distance, radii, level, max_level);
        if (w == 0.0){continue; }
        double force_common_part = r -> G * w / (r2 * distance);
        double fi = force_common_part * a->m;
        double fj = force_common_part * b->m;
        a->vx += dt * fj * dx;
        a->vy += dt * fj * dy;
        a->vz += dt * fj * dz;
        b->vx -= dt * fi * dx;
        b->vy -= dt * fi * dy;
        b->vz -= dt * fi * dz;
    }
}

static int need_deeper_level(const struct reb_particle* a,const struct reb_particle* b,double radius, double dt){
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double dz = b->z - a->z;
    double r2 = dx*dx + dy*dy + dz*dz;

    if (r2 <= radius * radius) {
        return 1;
    }
    double direction = 1.0;
    if (dt < 0)
    {
        direction = -1.0;
    }
    double dvx = (b->vx - a->vx) * direction;
    double dvy = (b->vy - a->vy) * direction;
    double dvz = (b->vz - a->vz) * direction;
    double rv = dx * dvx + dy * dvy + dz * dvz;
    double v2 = dvx * dvx + dvy * dvy + dvz * dvz;
    if(rv >= 0.0){
        return 0;
    }
    double t = fmin(-rv/v2,fabs(dt));
    const double x_near = dx + dvx*t;
    const double y_near = dy + dvy*t;
    const double z_near = dz + dvz*t;
    const double r2_min = x_near*x_near + y_near*y_near + z_near*z_near;
    return r2_min <= radius*radius;
}



static void evolve_level(struct reb_simulation* r,struct reb_particle* dh_particles, const int* particles, const struct SyMBA_Pair* pairs, size_t num_pairs, const double* radii,unsigned int level, unsigned int max_level, double dt){
    int* deeper_particles = NULL;
    struct SyMBA_Pair* deeper_pairs = NULL;
    size_t num_deeper = 0;
    level_wise_interaction_kick(r, dh_particles, pairs, num_pairs, radii, level, max_level, dt / 2);
    if(level < max_level){
        for(size_t k = 0; k < num_pairs; k++){
            const struct SyMBA_Pair pair = pairs[k];
            if (!need_deeper_level(&dh_particles[pair.i], &dh_particles[pair.j], radii[level], dt)) {
                continue;
            }
            if(num_deeper == 0){
                deeper_particles = calloc(r->N, sizeof(*deeper_particles));
                deeper_pairs = malloc(sizeof(struct SyMBA_Pair) * num_pairs);
            }
            deeper_particles[pair.i] = 1;
            deeper_particles[pair.j] = 1;
            deeper_pairs[num_deeper++] = pair;
        }
    }

    const double mu = r->G * r->particles[0].m;
    for(size_t i = 1; i < r->N; i++){
        if(particles[i] && (deeper_particles == NULL || !deeper_particles[i])){
            reb_integrator_whfast_kepler_solver(&dh_particles[i], mu, dt, r);

        }

    }

    if(num_deeper > 0){
        for (size_t i = 0; i < 3; i++) {
            evolve_level( r, dh_particles, deeper_particles, deeper_pairs, num_deeper, radii, level + 1, max_level, dt/3);
        }
    }


    free(deeper_particles);
    free(deeper_pairs);

    level_wise_interaction_kick(r, dh_particles, pairs, num_pairs, radii, level, max_level, dt / 2);
}


static void evolve(struct reb_simulation* r,struct reb_particle* dh_particles, const double* radii,unsigned int max_level, double dt){
    if (r->N <= 1) {
        return; // No particles to evolve
    }
    const size_t n = r->N - 1; // Exclude the central particle
    size_t a,b;

    if(n % 2 == 0){
        a = n / 2;
        b = n - 1;
    }else{ 
        a = n;
        b = (n - 1) / 2;
    }

    struct SyMBA_Pair* pairs;
    if(n==1){
        pairs = NULL;
    }else{
        pairs = malloc(sizeof(struct SyMBA_Pair) * a * b);
    }

    int* particles = calloc(r->N, sizeof(*particles));
    size_t k = 0;
    for(size_t i = 1; i < r->N ; i++){
        particles[i] = 1;
        for (size_t j = i + 1; j < r->N ; j++)
        {
            struct SyMBA_Pair pair = {i, j};
            pairs[k++] = pair;
        }

    }

    evolve_level(r, dh_particles, particles, pairs, k, radii, 0, max_level, dt);

    free(pairs);
    free(particles);
}


static void symba_step(struct reb_simulation* r, void* state){
    const struct reb_integrator_symba_state* s = state;

    if (!s || s->levels == 0 || !(s->r_outer > 0.0) || !isfinite(s->r_outer) ||
            (s->levels > 1 && (!(s->ratio > 1.0) || !isfinite(s->ratio)))){
        fail_step(r, "Invalid fixed-depth SyMBA shell configuration.");
        return;
    }
    const size_t n = r->N;
    if (!n || !(r->particles[0].m > 0.0) || !(r->G > 0.0) ||
            (r->N_active != SIZE_MAX && r->N_active != n) || r->N_var || r->additional_forces || r->gravity != REB_GRAVITY_BASIC){
        fail_step(r, "SyMBA requires active particles and isolated basic Newtonian gravity.");
        return;
    }

    if (r->dt == 0.0){
        r->dt_last_done = 0.0;
        return;}

    double* radii = malloc(sizeof(double) * s->levels);
    struct reb_particle* dh_particles = malloc(sizeof(struct reb_particle) * n);

    radii[0] = s->r_outer;
    for(unsigned int i = 1; i < s->levels; i++){
        radii[i] = radii[i - 1] / s->ratio;
        // very rare case, just to prevent underflow
        if (!(radii[i] > 0.0 && radii[i] < radii[i-1])){
            fail_step(r, "Shell radii must remain positive and strictly decreasing.");
            goto cleanup;
        }
    }

    const double dt = r->dt;

    reb_transformations_inertial_to_democraticheliocentric_posvel(r->particles, dh_particles, n, n);

    sun_drift(r, dh_particles, dt / 2);

    evolve(r, dh_particles, radii, s->levels - 1, dt);

    sun_drift(r, dh_particles, 0.5*dt);

    // Advance the barycenter, corresponding to the operation on Q[0] in Python.
    reb_integrator_whfast_com_step(r, dh_particles, dt);

    reb_transformations_democraticheliocentric_to_inertial_posvel(r->particles, dh_particles, n, n);

    r->t += dt;
    r->dt_last_done = dt;


cleanup:
    free(radii);
    free(dh_particles);

    return;
}

static const struct reb_binarydata_field_descriptor fields[] = {
    {
        "Number of levels including level zero.", REB_UINT,
        "levels", offsetof(struct reb_integrator_symba_state, levels), 0, 0, 0
    },
    {
        "Outermost radius in simulation length units.", REB_DOUBLE,
        "r_outer", offsetof(struct reb_integrator_symba_state, r_outer), 0, 0, 0
    },
    {
        "Ratio of successive radii.", REB_DOUBLE,
        "ratio", offsetof(struct reb_integrator_symba_state, ratio), 0, 0, 0
    },
    {0},
};


const struct reb_integrator reb_integrator_symba = {
    .documentation = "SyMBA with geometric shells and three substeps per level.",
    .step = symba_step,
    .create = create_symba_state,
    .free = free,
    .field_descriptor_list = fields,
};
