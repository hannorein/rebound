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
    if (whfast == NULL){
        return;
    }
    free(whfast->nodes);
    free(whfast->barycenters);
    free(whfast->masses);
    free(whfast->tree);
    free(whfast);
}

// Resolve a subtree to its shared barycenter. A massless attachment aliases
// its massive child, so a chain of test particles adds no barycenter work.
static size_t reb_integrator_whfast_hj_barycenter_index( const struct reb_integrator_whfast_hj_state* const whfast, const size_t index){
    return index < whfast->tree_N ? index : whfast->nodes[index - whfast->tree_N].barycenter;
}

static struct reb_particle* reb_integrator_whfast_hj_barycenter( struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const size_t index){
    const size_t host = reb_integrator_whfast_hj_barycenter_index(whfast, index);
    return host < whfast->tree_N ? &r->particles[host] : &whfast->barycenters[host - whfast->tree_N];
}

static double reb_integrator_whfast_hj_subtree_mass( const struct reb_integrator_whfast_hj_state* const whfast, const size_t index){
    return index < whfast->tree_N ? whfast->masses[index] : whfast->nodes[index - whfast->tree_N].q.m;
}

static void reb_integrator_whfast_hj_initialize_internal( struct reb_integrator_whfast_hj_state* const whfast, const size_t index, const size_t primary, const size_t secondary){
    struct reb_integrator_whfast_hj_node* const node = &whfast->nodes[index - whfast->tree_N];
    const double primary_mass = reb_integrator_whfast_hj_subtree_mass(whfast, primary);
    const double secondary_mass = reb_integrator_whfast_hj_subtree_mass(whfast, secondary);
    const double total_mass = primary_mass + secondary_mass;
    node->primary = primary;
    node->secondary = secondary;
    node->primary_offset = secondary_mass/total_mass;
    node->secondary_offset = primary_mass/total_mass;
    node->q.m = total_mass;
    if (primary_mass == 0.){
        node->barycenter = reb_integrator_whfast_hj_barycenter_index(whfast, secondary);
    }else if (secondary_mass == 0.){
        node->barycenter = reb_integrator_whfast_hj_barycenter_index(whfast, primary);
    }else{
        node->barycenter = whfast->tree_N + whfast->N_barycenters++;
    }
}

struct reb_integrator_whfast_hj_tree_parser {
    struct reb_integrator_whfast_hj_state* whfast;
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
// The caller supplies valid syntax, each particle exactly once, and positive
// total mass for every binary. Whitespace is allowed between tokens.
static size_t reb_integrator_whfast_hj_parse_tree_node(struct reb_integrator_whfast_hj_tree_parser* const parser) {
    reb_integrator_whfast_hj_tree_parser_skip_space(parser);

    if (*parser->cursor == '['){
        parser->cursor++;
        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        if (*parser->cursor == ']'){
            parser->cursor++;
            return SIZE_MAX;
        }

        size_t children[2];
        for (int i=0; i<2; i++){
            children[i] = reb_integrator_whfast_hj_parse_tree_node(parser);
            reb_integrator_whfast_hj_tree_parser_skip_space(parser);
            parser->cursor++;  // Consume ',' after the first child, ']' after the second.
        }

        // Unique in-range leaves allow at most tree_N - 1 internal binary nodes.
        const size_t index = parser->next_internal++;
        // here
        reb_integrator_whfast_hj_initialize_internal(parser->whfast, index, children[0], children[1]);
        return index;
    }

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
    const size_t orbit_count = r->N - 1;
    whfast->masses = calloc(r->N, sizeof(*whfast->masses));
    whfast->nodes = calloc(orbit_count, sizeof(*whfast->nodes));
    // At most K-1 new barycenters are needed for K nonzero-mass particles.
    size_t massive_count = 0;
    for (size_t i=0; i<r->N; i++){
        massive_count += r->particles[i].m != 0.;
        whfast->masses[i] = r->particles[i].m;
    }
    whfast->barycenters = calloc(massive_count - 1, sizeof(*whfast->barycenters));

    struct reb_integrator_whfast_hj_tree_parser parser = {
        .whfast = whfast,
        .cursor = tree,
        .next_internal = r->N,
    };
    reb_integrator_whfast_hj_parse_tree_node(&parser);

    return 0;
}

// Positions/velocities and accelerations have separate passes: the gravity
// evaluation changes only accelerations, so there is no midpoint round trip
// for the already-current Jacobi positions and velocities.
static void reb_integrator_whfast_hj_from_inertial( struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const int accelerations){
    for (size_t i=0; i+1<whfast->tree_N; i++){
        struct reb_integrator_whfast_hj_node* const node = &whfast->nodes[i];
        const struct reb_particle* const primary = reb_integrator_whfast_hj_barycenter(r, whfast, node->primary);
        const struct reb_particle* const secondary = reb_integrator_whfast_hj_barycenter(r, whfast, node->secondary);
        struct reb_particle* const q = &node->q;
        struct reb_particle* const com = reb_integrator_whfast_hj_barycenter(r, whfast, whfast->tree_N + i);
        const int shared = com == primary || com == secondary;
#define HJ_FORWARD(field) \
        q->field = secondary->field - primary->field; \
        if (!shared) com->field = node->secondary_offset*primary->field + node->primary_offset*secondary->field;
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
        struct reb_particle* const primary = reb_integrator_whfast_hj_barycenter(r, whfast, node->primary);
        struct reb_particle* const secondary = reb_integrator_whfast_hj_barycenter(r, whfast, node->secondary);
        const struct reb_particle* const com = reb_integrator_whfast_hj_barycenter(r, whfast, whfast->tree_N + i - 1);
        const struct reb_particle* const q = &node->q;
        // Shared hosts have already been reconstructed by their ancestors.
        // Only place the massless child; do not rewrite the host for each tracer.
#define HJ_PRIMARY(field) primary->field = com->field - node->primary_offset*q->field
#define HJ_SECONDARY(field) secondary->field = com->field + node->secondary_offset*q->field
        if (primary != com){
            HJ_PRIMARY(x); HJ_PRIMARY(y); HJ_PRIMARY(z);
            HJ_PRIMARY(vx); HJ_PRIMARY(vy); HJ_PRIMARY(vz);
        }
        if (secondary != com){
            HJ_SECONDARY(x); HJ_SECONDARY(y); HJ_SECONDARY(z);
            HJ_SECONDARY(vx); HJ_SECONDARY(vy); HJ_SECONDARY(vz);
        }
#undef HJ_PRIMARY
#undef HJ_SECONDARY
    }
}

static void reb_integrator_whfast_hj_interaction_step( const struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const double dt){
    for (size_t i=0; i+1<whfast->tree_N; i++){
        struct reb_particle* const p = &whfast->nodes[i].q;
        p->vx += dt*p->ax;
        p->vy += dt*p->ay;
        p->vz += dt*p->az;

        const double inverse_r2 = 1./(p->x*p->x + p->y*p->y + p->z*p->z);
        const double inverse_r = sqrt(inverse_r2);
        const double prefactor = dt*r->G*p->m*inverse_r*inverse_r2;
        p->vx += prefactor*p->x;
        p->vy += prefactor*p->y;
        p->vz += prefactor*p->z;
    }
}

static void reb_integrator_whfast_hj_kepler_step( const struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const double dt){
    for (size_t i=0; i+1<whfast->tree_N; i++){
        struct reb_particle* const q = &whfast->nodes[i].q;
        reb_integrator_whfast_kepler_solver(q, q->m*r->G, dt, r);
    }
}

static void reb_integrator_whfast_hj_com_step( struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast, const double dt){
    const size_t node_count = 2*whfast->tree_N-1;
    struct reb_particle* const com = reb_integrator_whfast_hj_barycenter(r, whfast, node_count - 1);
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
    reb_integrator_whfast_hj_com_step(r, whfast, dt/2.);
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
    reb_integrator_whfast_hj_com_step(r, whfast, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast);

    r->t += dt;
    r->dt_last_done = dt;
}
