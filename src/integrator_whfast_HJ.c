#include "rebound.h"
#include "rebound_internal.h"
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "integrator_whfast.h"
#include "integrator_whfast_HJ.h"
#include "binarydata.h"

void* reb_integrator_whfast_hj_create(void);
void reb_integrator_whfast_hj_free(void* state);
void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state);
const struct reb_binarydata_field_descriptor reb_integrator_whfast_hj_field_descriptor_list[];

const struct reb_integrator reb_integrator_whfast_hj = {
    .documentation =
    "WHFast HJ is a hierarchical Jacobi-coordinate variant of WHFast. "
    "It requires a fixed user-supplied binary hierarchy, which is compiled once "
    "into compact orbit arrays with shared barycenters for massless particles.",
    .step = reb_integrator_whfast_hj_step,
    .create = reb_integrator_whfast_hj_create,
    .free = reb_integrator_whfast_hj_free,
    .field_descriptor_list = reb_integrator_whfast_hj_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_whfast_hj_field_descriptor_list[] = {
    { "Tree structure. Format example: [[1,2],3]. The tree can only be set once. "
        "To change the tree aftwards, you need to first reset the integrator.",
        REB_STRING, "tree", offsetof(struct reb_integrator_whfast_hj_state, tree), 0, 0, 0},
    { 0 },
};

void* reb_integrator_whfast_hj_create(void) {
    return calloc(1, sizeof(struct reb_integrator_whfast_hj_state));
}

void reb_integrator_whfast_hj_free(void* state) {
    struct reb_integrator_whfast_hj_state* const whfast = state;
    free(whfast->nodes);
    free(whfast->barycenters);
    free(whfast->tree);
    free(whfast);
}

static struct reb_particle* reb_integrator_whfast_hj_get_particle( struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const size_t index){
    if (index < whfast->tree_N){
        // Real particle
        return &r->particles[index];
    }else{
        // Barycenter particle
        size_t b_index = whfast->nodes[index - whfast->tree_N].barycenter;
        return &whfast->barycenters[b_index];
    }
}

struct reb_integrator_whfast_hj_tree_parser {
    const char* cursor;
    size_t next_internal;
};

static void reb_integrator_whfast_hj_tree_parser_skip_space(struct reb_integrator_whfast_hj_tree_parser* const parser) {
    while (isspace((unsigned char)*parser->cursor)){
        parser->cursor++;
    }
}

// The input grammar is recursive, but the result is a flat postordered array.
// Recursion is used only while compiling a user-supplied tree, never per timestep.
static size_t reb_integrator_whfast_hj_parse_tree_node(const struct reb_simulation* r, struct reb_integrator_whfast_hj_state* const whfast, struct reb_integrator_whfast_hj_tree_parser* const parser) {
    reb_integrator_whfast_hj_tree_parser_skip_space(parser);

    if (*parser->cursor == '['){
        parser->cursor++;
        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        size_t primary   = reb_integrator_whfast_hj_parse_tree_node(r, whfast, parser);
        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        parser->cursor++;  // Consume ','
        size_t secondary = reb_integrator_whfast_hj_parse_tree_node(r, whfast, parser);
        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        parser->cursor++;  // Consume ']'

        // Unique in-range leaves allow at most tree_N - 1 internal binary nodes.
        const size_t index = parser->next_internal++;
        struct reb_integrator_whfast_hj_node* const node = &whfast->nodes[index - whfast->tree_N];
        node->primary = primary;
        node->secondary = secondary;
        node->barycenter = whfast->N_barycenters++;
        return index;
    }

    // Node is a real particle
    size_t particle_number = 0;
    while (isdigit((unsigned char)*parser->cursor)){
        particle_number = 10*particle_number + (size_t)(*parser->cursor - '0');
        parser->cursor++;
    }
    return particle_number - 1;
}

int reb_integrator_whfast_hj_set_tree(struct reb_simulation* const r, const char* const tree) {
    struct reb_integrator_whfast_hj_state* const whfast = r->integrator.state;

    whfast->tree_N = r->N;
    whfast->nodes = calloc(r->N-1, sizeof(*whfast->nodes));
    whfast->barycenters = calloc(r->N-1, sizeof(*whfast->barycenters));

    struct reb_integrator_whfast_hj_tree_parser parser = {
        .cursor = tree,
        .next_internal = r->N,
    };
    reb_integrator_whfast_hj_parse_tree_node(r, whfast, &parser);

    return 0;
}

// Positions/velocities and accelerations have separate passes: the gravity
// evaluation changes only accelerations, so there is no midpoint round trip
// for the already-current Jacobi positions and velocities.
static void reb_integrator_whfast_hj_from_inertial( struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const int accelerations){
    for (size_t i=0; i+1<whfast->tree_N; i++){
        struct reb_integrator_whfast_hj_node* const node = &whfast->nodes[i];
        const struct reb_particle* const primary = reb_integrator_whfast_hj_get_particle(r, whfast, node->primary);
        const struct reb_particle* const secondary = reb_integrator_whfast_hj_get_particle(r, whfast, node->secondary);
        struct reb_particle* const q = &node->q;
        struct reb_particle* const barycenter = &whfast->barycenters[node->barycenter];
        barycenter->m = primary->m + secondary->m;
#define HJ_FORWARD(field) \
        q->field = secondary->field - primary->field; \
        barycenter->field = (primary->m*primary->field + secondary->m*secondary->field)/barycenter->m;
        if (accelerations){
            HJ_FORWARD(ax); HJ_FORWARD(ay); HJ_FORWARD(az);
        }else{
            HJ_FORWARD(x); HJ_FORWARD(y); HJ_FORWARD(z);
            HJ_FORWARD(vx); HJ_FORWARD(vy); HJ_FORWARD(vz);
        }
#undef HJ_FORWARD
    }
}

static void reb_integrator_whfast_hj_to_inertial( struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast){
    for (size_t i=whfast->tree_N-1; i>0;i--){
        const struct reb_integrator_whfast_hj_node* const node = &whfast->nodes[i-1];
        struct reb_particle* const primary = reb_integrator_whfast_hj_get_particle(r, whfast, node->primary);
        struct reb_particle* const secondary = reb_integrator_whfast_hj_get_particle(r, whfast, node->secondary);
        const struct reb_particle* const barycenter = &whfast->barycenters[node->barycenter];
        const struct reb_particle* const q = &node->q;
#define HJ_PRIMARY(field) primary->field = barycenter->field - secondary->m/barycenter->m*q->field
#define HJ_SECONDARY(field) secondary->field = barycenter->field + primary->m/barycenter->m*q->field
        HJ_PRIMARY(x); HJ_PRIMARY(y); HJ_PRIMARY(z);
        HJ_PRIMARY(vx); HJ_PRIMARY(vy); HJ_PRIMARY(vz);
        HJ_SECONDARY(x); HJ_SECONDARY(y); HJ_SECONDARY(z);
        HJ_SECONDARY(vx); HJ_SECONDARY(vy); HJ_SECONDARY(vz);
#undef HJ_PRIMARY
#undef HJ_SECONDARY
    }
}

static void reb_integrator_whfast_hj_interaction_step( const struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const double dt){
    for (size_t i=0; i+1<whfast->tree_N; i++){
        struct reb_integrator_whfast_hj_node* node = &whfast->nodes[i];
        struct reb_particle* const q = &node->q;
        struct reb_particle* const barycenter = &whfast->barycenters[node->barycenter];
        q->vx += dt*q->ax;
        q->vy += dt*q->ay;
        q->vz += dt*q->az;

        const double inverse_r2 = 1./(q->x*q->x + q->y*q->y + q->z*q->z);
        const double inverse_r = sqrt(inverse_r2);
        const double prefactor = dt*r->G*barycenter->m*inverse_r*inverse_r2;
        q->vx += prefactor*q->x;
        q->vy += prefactor*q->y;
        q->vz += prefactor*q->z;
    }
}

static void reb_integrator_whfast_hj_kepler_step( const struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const double dt){
    for (size_t i=0; i+1<whfast->tree_N; i++){
        struct reb_integrator_whfast_hj_node* node = &whfast->nodes[i];
        struct reb_particle* barycenter = &whfast->barycenters[node->barycenter];
        reb_integrator_whfast_kepler_solver(&node->q, barycenter->m*r->G, dt, r);
    }
}

static void reb_integrator_whfast_hj_com_step(struct reb_integrator_whfast_hj_state* const whfast, const double dt){
    struct reb_particle* const com = &whfast->barycenters[whfast->tree_N-2];
    com->x += dt*com->vx;
    com->y += dt*com->vy;
    com->z += dt*com->vz;
}

void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state) {
    struct reb_integrator_whfast_hj_state* const whfast = state;
    const double dt = r->dt;

    if (!whfast->tree_N){
        if (!whfast->tree){
            reb_simulation_error(r, "WHFast HJ requires a fixed tree to be set.");
            return;
        }
        reb_integrator_whfast_hj_set_tree(r, whfast->tree);
    }

    if (whfast->tree_N != r->N){
        reb_simulation_error(r, "WHFast HJ particle count changed after the fixed tree was set.");
        return;
    }
    reb_integrator_whfast_hj_from_inertial(r, whfast, 0);

    reb_integrator_whfast_hj_kepler_step(r, whfast, dt/2.);
    reb_integrator_whfast_hj_com_step(whfast, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast);

    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    const int refresh_coordinates = r->additional_forces || r->gravity == REB_GRAVITY_CUSTOM;
    reb_simulation_update_acceleration(r);
    // Preserve the old behavior even for callbacks that also modify x/v.
    if (refresh_coordinates){
        reb_integrator_whfast_hj_from_inertial(r, whfast, 0);
    }
    reb_integrator_whfast_hj_from_inertial(r, whfast, 1);
    reb_integrator_whfast_hj_interaction_step(r, whfast, dt);

    reb_integrator_whfast_hj_kepler_step(r, whfast, dt/2.);
    reb_integrator_whfast_hj_com_step(whfast, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast);

    r->t += dt;
    r->dt_last_done = dt;
}
