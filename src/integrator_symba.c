#include "rebound.h"
#include "integrator_whfast.h"
#include "transformations.h"
#include "integrator_symba.h"
#include "assert.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

typedef struct reb_integrator_symba_state SymbaState;

static void* create_symba_state(void) {
    SymbaState *state = malloc(sizeof(SymbaState));
    if(state != NULL) {
        state->levels = 3;
        state->r_outer = 0.15;
        state->ratio = 2.08;
    }
    return state;

    
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

static bool level_wise_interaction_kick(const struct reb_simulation* r, struct reb_particle* dh_particles,
    const double* radii,unsigned int level, unsigned int max_level,double dt){
        for(size_t i = 1; i < r->N; i++){
            for(size_t j = i+1; j < r->N; j++){
                double dx = dh_particles[j].x - dh_particles[i].x;
                double dy = dh_particles[j].y - dh_particles[i].y;
                double dz = dh_particles[j].z - dh_particles[i].z;
                double r2 = dx*dx + dy*dy + dz*dz;
                if(r2 == 0.0){
                    return false;
                }
                const double distance = sqrt(r2);
                const double w = level_weight(distance, radii, level, max_level);
                double force_common_part = r -> G * w / (r2 * distance);
                double fi = force_common_part * dh_particles[i].m;
                double fj = force_common_part * dh_particles[j].m;
                dh_particles[i].vx += dt * fj * dx;
                dh_particles[i].vy += dt * fj * dy;
                dh_particles[i].vz += dt * fj * dz;
                dh_particles[j].vx -= dt * fi * dx;
                dh_particles[j].vy -= dt * fi * dy;
                dh_particles[j].vz -= dt * fi * dz;
            }
        }
    return true;
}


static bool evolve_level(struct reb_simulation* r, struct reb_particle* dh_particles,
    const double* radii,unsigned int level, unsigned int max_level,double dt){
        if (!level_wise_interaction_kick(r,dh_particles,radii,level,max_level,dt/2)){
            return false;
        }
        if(level == max_level){
            reb_integrator_whfast_kepler_step(r, dh_particles,REB_INTEGRATOR_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC,dt);
        }else{
            for(unsigned int i = 0; i < 3; i++){
                if(!evolve_level(r, dh_particles, radii, level + 1, max_level, dt / 3)){
                    return false;
                }
            }
        }
        return level_wise_interaction_kick(r,dh_particles,radii,level,max_level,dt/2);
    }
    
    
static void fail_step(struct reb_simulation* r, const char* message){
    r->status = REB_STATUS_GENERIC_ERROR;
    reb_simulation_error(r, message);
}

static void symba_step(struct reb_simulation* r, void* state){
    const SymbaState* s = state;

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

    double* radii = malloc(sizeof(double) * s->levels);
    struct reb_particle* dh_particles = malloc(sizeof(struct reb_particle) * n);

    if(!radii || !dh_particles){
        fail_step(r, "Memory allocation failed.");
        goto cleanup;
    }


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

    if(!evolve_level(r, dh_particles, radii, 0, s->levels - 1, dt)){
        fail_step(r, "Integration failed due to a  collision.");
        goto cleanup;
    }

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
        "levels", offsetof(SymbaState, levels), 0, 0, 0
    },
    {
        "Outermost radius in simulation length units.", REB_DOUBLE,
        "r_outer", offsetof(SymbaState, r_outer), 0, 0, 0
    },
    {
        "Ratio of successive radii.", REB_DOUBLE,
        "ratio", offsetof(SymbaState, ratio), 0, 0, 0
    },
    {0},
};
              

const struct reb_integrator reb_integrator_symba = {
    .documentation = "Fixed-depth SyMBA with geometric shells and three substeps per level.",
    .step = symba_step,
    .create = create_symba_state,
    .free = free,
    .field_descriptor_list = fields,
};