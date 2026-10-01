
#include "rebound.h"
#include "rebound_internal.h"
#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include "transformations.h"
#include "particle.h"
#include "tools.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator_whfast_HJ.h"
#include "binarydata.h"

#define WHFAST_HJ_UNBOUND_DISTANCE_FACTOR 0.3

void* reb_integrator_whfast_hj_create();
void reb_integrator_whfast_hj_free(void* state);
void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state);
static void reb_integrator_whfast_hj_node_free(struct hj_node* const node);
const struct reb_binarydata_field_descriptor reb_integrator_whfast_hj_field_descriptor_list[];

const struct reb_integrator reb_integrator_whfast_hj = {
    .documentation =
    "WHFast HJ is a hierarchical Jacobi-coordinate variant of WHFast. "
    "By default it rebuilds the binary hierarchy tree at every timestep. "
    "A user-supplied tree can be cached with sim.integrate(..., given_tree=True, tree=...).",
    .step = reb_integrator_whfast_hj_step,
    .create = reb_integrator_whfast_hj_create,
    .free = reb_integrator_whfast_hj_free,
    .field_descriptor_list = reb_integrator_whfast_hj_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_whfast_hj_field_descriptor_list[] = {
    { "Use a user-supplied HJ tree instead of rebuilding the tree each timestep.",
        REB_INT,        "given_tree",          offsetof(struct reb_integrator_whfast_hj_state, given_tree), 0, 0, 0},
    { 0 }, // Null terminated list
};

// create and free
void* reb_integrator_whfast_hj_create(){
    // Allocate memory and set default parameters.
    struct reb_integrator_whfast_hj_state* whfast = calloc(sizeof(struct reb_integrator_whfast_hj_state),1);
    return whfast;
}

static void reb_integrator_whfast_hj_node_free(struct hj_node *const node)
{
    if (node == NULL){
        return;
    }

    size_t stack_size = 1;
    size_t stack_capacity = 128;
    struct hj_node** stack = malloc(stack_capacity*sizeof(struct hj_node*));
    if (stack == NULL){
        return;
    }
    stack[0] = node;

    while (stack_size > 0){
        struct hj_node* const current = stack[--stack_size];
        if (current->primary != NULL){
            if (stack_size == stack_capacity){
                stack_capacity *= 2;
                struct hj_node** const new_stack = realloc(stack, stack_capacity*sizeof(struct hj_node*));
                if (new_stack == NULL){
                    free(stack);
                    return;
                }
                stack = new_stack;
            }
            stack[stack_size++] = current->primary;
        }
        if (current->secondary != NULL){
            if (stack_size == stack_capacity){
                stack_capacity *= 2;
                struct hj_node** const new_stack = realloc(stack, stack_capacity*sizeof(struct hj_node*));
                if (new_stack == NULL){
                    free(stack);
                    return;
                }
                stack = new_stack;
            }
            stack[stack_size++] = current->secondary;
        }
        free(current);
    }
    free(stack);
}

void reb_integrator_whfast_hj_free(void* p){
    struct reb_integrator_whfast_hj_state* whfast = p;
    reb_integrator_whfast_hj_node_free(whfast->root);
    free(whfast);
}

static struct hj_node* reb_integrator_whfast_hj_node_create_leaf(const struct reb_particle particle, const int particle_index)
{
    struct hj_node* const node = calloc(1, sizeof(struct hj_node));
    if (node == NULL){
        return NULL;
    }
    node->barycenter_particle = particle;
    node->particle_index = particle_index;
    return node;
}

static struct hj_node* reb_integrator_whfast_hj_node_create_binary_ordered(struct hj_node* const primary, struct hj_node* const secondary)
{
    struct hj_node* const node = calloc(1, sizeof(struct hj_node));
    if (node == NULL){
        return NULL;
    }

    node->primary = primary;
    node->secondary = secondary;
    node->barycenter_particle = reb_particle_com_of_pair(node->primary->barycenter_particle, node->secondary->barycenter_particle);
    node->particle_index = -1;
    return node;
}

static struct hj_node* reb_integrator_whfast_hj_node_create_binary(struct hj_node* const a, struct hj_node* const b)
{
    if (a->barycenter_particle.m >= b->barycenter_particle.m){
        return reb_integrator_whfast_hj_node_create_binary_ordered(a, b);
    }
    return reb_integrator_whfast_hj_node_create_binary_ordered(b, a);
}

static void reb_integrator_whfast_hj_node_list_free(struct hj_node** const nodes, const size_t N)
{
    for (size_t i=0; i<N; i++){
        reb_integrator_whfast_hj_node_free(nodes[i]);
    }
    free(nodes);
}

static struct reb_orbit reb_integrator_whfast_hj_orbit_for_pair(const double G, const struct hj_node* const a, const struct hj_node* const b)
{
    if (a->barycenter_particle.m >= b->barycenter_particle.m){
        return reb_orbit_from_particle(G, b->barycenter_particle, a->barycenter_particle);
    }
    return reb_orbit_from_particle(G, a->barycenter_particle, b->barycenter_particle);
}

static double reb_integrator_whfast_hj_bound_pair_timescale(const double G, const struct hj_node* const a, const struct hj_node* const b)
{
    const struct reb_orbit o = reb_integrator_whfast_hj_orbit_for_pair(G, a, b);
    if (isfinite(o.P) && o.P > 0.){
        return o.P;
    }
    return INFINITY;
}

static double reb_integrator_whfast_hj_unbound_pair_timescale(const double G, const struct hj_node* const a, const struct hj_node* const b)
{
    const struct reb_orbit o = reb_integrator_whfast_hj_orbit_for_pair(G, a, b);
    const double M = a->barycenter_particle.m + b->barycenter_particle.m;
    const double mu = G*M;
    if (!isfinite(mu) || mu <= 0. || !isfinite(o.e) || o.e <= 1.){
        return INFINITY;
    }

    const struct reb_particle pa = a->barycenter_particle;
    const struct reb_particle pb = b->barycenter_particle;

    const double rx = pb.x - pa.x;
    const double ry = pb.y - pa.y;
    const double rz = pb.z - pa.z;
    const double vx = pb.vx - pa.vx;
    const double vy = pb.vy - pa.vy;
    const double vz = pb.vz - pa.vz;

    const double hx = ry*vz - rz*vy;
    const double hy = rz*vx - rx*vz;
    const double hz = rx*vy - ry*vx;
    const double h2 = hx*hx + hy*hy + hz*hz;
    const double denominator = mu*(1.0 + o.e);
    if (!isfinite(h2) || h2 <= 0. || !isfinite(denominator) || denominator <= 0.){
        return INFINITY;
    }

    const double rp = h2/denominator;
    if (!isfinite(rp) || rp <= 0.){
        return INFINITY;
    }

    const double tperi = sqrt(rp*rp*rp/denominator);
    return (isfinite(tperi) && tperi > 0.) ? tperi : INFINITY;
}

static double reb_integrator_whfast_hj_pair_distance(const struct hj_node* const a, const struct hj_node* const b)
{
    const double dx = b->barycenter_particle.x - a->barycenter_particle.x;
    const double dy = b->barycenter_particle.y - a->barycenter_particle.y;
    const double dz = b->barycenter_particle.z - a->barycenter_particle.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

//----------------------------------------------------------------------------
// build binary orbit hierarchy tree
static int reb_integrator_whfast_hj_build_tree(struct reb_simulation* const r, struct reb_integrator_whfast_hj_state* const whfast)
{
    reb_integrator_whfast_hj_node_free(whfast->root);
    whfast->root = NULL;
    whfast->tree_N = 0;

    size_t N = r->N;
    if (N == 0){
        return 0;
    }

    struct hj_node** const nodes = calloc(N, sizeof(struct hj_node*));
    if (nodes == NULL){
        return 1;
    }

    for (size_t i=0; i<N; i++){
        nodes[i] = reb_integrator_whfast_hj_node_create_leaf(r->particles[i], (int)i);
        if (nodes[i] == NULL){
            reb_integrator_whfast_hj_node_list_free(nodes, i);
            return 1;
        }
    }

    while (N > 1){
        size_t imin = SIZE_MAX;
        size_t jmin = SIZE_MAX;
        double Pmin = INFINITY;
        double dmin_bound = INFINITY;

        for (size_t i=1; i<N; i++){
            for (size_t j=0; j<i; j++){
                double P = reb_integrator_whfast_hj_bound_pair_timescale(r->G, nodes[i], nodes[j]);
                // if (P == INFINITY)
                // {
                //     double Pperi = reb_integrator_whfast_hj_unbound_pair_timescale(r->G, nodes[i], nodes[j]);
                //     if (Pperi < Pmin && reb_integrator_whfast_hj_pair_distance(nodes[i], nodes[j]) < WHFAST_HJ_UNBOUND_DISTANCE_FACTOR * dmin_bound)
                //     {
                //         Pmin = Pperi;
                //         imin = i;
                //         jmin = j;
                //     }
                // }
                if (P < Pmin){
                    Pmin = P;
                    imin = i;
                    jmin = j;
                    dmin_bound = reb_integrator_whfast_hj_pair_distance(nodes[i], nodes[j]);
                }
            }
        }

        for (size_t i=1; i<N; i++){
            for (size_t j=0; j<i; j++){
                const double Pperi = reb_integrator_whfast_hj_unbound_pair_timescale(r->G, nodes[i], nodes[j]);
                if (Pperi < Pmin && reb_integrator_whfast_hj_pair_distance(nodes[i], nodes[j]) < WHFAST_HJ_UNBOUND_DISTANCE_FACTOR*dmin_bound){
                    Pmin = Pperi;
                    imin = i;
                    jmin = j;
                }
            }
        }

        if (Pmin == INFINITY){
            reb_integrator_whfast_hj_node_list_free(nodes, N);
            return 1;
        }

        struct hj_node* const merged = reb_integrator_whfast_hj_node_create_binary(nodes[imin], nodes[jmin]);
        if (merged == NULL){
            reb_integrator_whfast_hj_node_list_free(nodes, N);
            return 1;
        }

        nodes[imin] = merged;
        nodes[jmin] = nodes[N-1];
        N--;
    }

    whfast->root = nodes[0];
    whfast->tree_N = r->N;
    free(nodes);
    return 0;
}

//----------------------------------------------------------------------------
// Tree-Jacobi coordinate transforms.
static int reb_integrator_whfast_hj_node_is_leaf(const struct hj_node* const node)
{
    return node != NULL && node->primary == NULL && node->secondary == NULL;
}

static int reb_integrator_whfast_hj_node_stack_push(struct hj_node*** stack, size_t* const stack_size, size_t* const stack_capacity, struct hj_node* const node)
{
    if (node == NULL){
        return 0;
    }
    if (*stack_size == *stack_capacity){
        *stack_capacity *= 2;
        struct hj_node** const new_stack = realloc(*stack, *stack_capacity*sizeof(struct hj_node*));
        if (new_stack == NULL){
            return 1;
        }
        *stack = new_stack;
    }
    (*stack)[(*stack_size)++] = node;
    return 0;
}

struct reb_integrator_whfast_hj_visit_stack_item {
    struct hj_node* node;
    int visited;
};

static int reb_integrator_whfast_hj_visit_stack_push(
    struct reb_integrator_whfast_hj_visit_stack_item** stack,
    size_t* const stack_size,
    size_t* const stack_capacity,
    struct hj_node* const node,
    const int visited
){
    if (node == NULL){
        return 0;
    }
    if (*stack_size == *stack_capacity){
        *stack_capacity *= 2;
        struct reb_integrator_whfast_hj_visit_stack_item* const new_stack = realloc(*stack, *stack_capacity*sizeof(struct reb_integrator_whfast_hj_visit_stack_item));
        if (new_stack == NULL){
            return 1;
        }
        *stack = new_stack;
    }
    (*stack)[(*stack_size)++] = (struct reb_integrator_whfast_hj_visit_stack_item){node, visited};
    return 0;
}

static void reb_integrator_whfast_hj_tree_append(char* const buffer, const size_t buffer_size, size_t* const required, const char* const text)
{
    const size_t len = strlen(text);
    if (buffer != NULL && buffer_size > 0 && *required < buffer_size - 1){
        size_t copy_len = buffer_size - 1 - *required;
        if (copy_len > len){
            copy_len = len;
        }
        memcpy(buffer + *required, text, copy_len);
        buffer[*required + copy_len] = '\0';
    }
    *required += len;
}

static void reb_integrator_whfast_hj_tree_node_to_string(const struct hj_node* const node, char* const buffer, const size_t buffer_size, size_t* const required)
{
    if (reb_integrator_whfast_hj_node_is_leaf(node)){
        char leaf[32];
        snprintf(leaf, sizeof(leaf), "%d", node->particle_index + 1);
        reb_integrator_whfast_hj_tree_append(buffer, buffer_size, required, leaf);
        return;
    }

    reb_integrator_whfast_hj_tree_append(buffer, buffer_size, required, "[");
    reb_integrator_whfast_hj_tree_node_to_string(node->primary, buffer, buffer_size, required);
    reb_integrator_whfast_hj_tree_append(buffer, buffer_size, required, ",");
    reb_integrator_whfast_hj_tree_node_to_string(node->secondary, buffer, buffer_size, required);
    reb_integrator_whfast_hj_tree_append(buffer, buffer_size, required, "]");
}

REB_API int reb_integrator_whfast_hj_tree_to_string(struct reb_simulation* const r, char* const buffer, const size_t buffer_size)
{
    if (buffer != NULL && buffer_size > 0){
        buffer[0] = '\0';
    }
    if (r == NULL){
        return -1;
    }

    struct reb_integrator_whfast_hj_state whfast = {0};
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        return -1;
    }

    size_t required = 0;
    if (whfast.root == NULL){
        reb_integrator_whfast_hj_tree_append(buffer, buffer_size, &required, "[]");
    }else{
        reb_integrator_whfast_hj_tree_node_to_string(whfast.root, buffer, buffer_size, &required);
    }

    reb_integrator_whfast_hj_node_free(whfast.root);
    if (required > (size_t)INT_MAX){
        return -1;
    }
    return (int)required;
}

struct reb_integrator_whfast_hj_tree_parser {
    struct reb_simulation* r;
    const char* cursor;
    unsigned char* used;
    size_t leaf_count;
    const char* error;
};

static void reb_integrator_whfast_hj_tree_parser_skip_space(struct reb_integrator_whfast_hj_tree_parser* const parser)
{
    while (isspace((unsigned char)*parser->cursor)){
        parser->cursor++;
    }
}

static void reb_integrator_whfast_hj_tree_parser_error(struct reb_integrator_whfast_hj_tree_parser* const parser, const char* const error)
{
    if (parser->error == NULL){
        parser->error = error;
    }
}

static struct hj_node* reb_integrator_whfast_hj_parse_tree_node(struct reb_integrator_whfast_hj_tree_parser* const parser)
{
    reb_integrator_whfast_hj_tree_parser_skip_space(parser);

    if (*parser->cursor == '['){
        parser->cursor++;
        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        if (*parser->cursor == ']'){
            parser->cursor++;
            return NULL;
        }

        struct hj_node* const primary = reb_integrator_whfast_hj_parse_tree_node(parser);
        if (parser->error != NULL){
            return NULL;
        }
        if (primary == NULL){
            reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: empty subtree.");
            return NULL;
        }

        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        if (*parser->cursor != ','){
            reb_integrator_whfast_hj_node_free(primary);
            reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: expected ','.");
            return NULL;
        }
        parser->cursor++;

        struct hj_node* const secondary = reb_integrator_whfast_hj_parse_tree_node(parser);
        if (parser->error != NULL){
            reb_integrator_whfast_hj_node_free(primary);
            return NULL;
        }
        if (secondary == NULL){
            reb_integrator_whfast_hj_node_free(primary);
            reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: empty subtree.");
            return NULL;
        }

        reb_integrator_whfast_hj_tree_parser_skip_space(parser);
        if (*parser->cursor != ']'){
            reb_integrator_whfast_hj_node_free(primary);
            reb_integrator_whfast_hj_node_free(secondary);
            reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: expected ']'.");
            return NULL;
        }
        parser->cursor++;

        struct hj_node* const node = reb_integrator_whfast_hj_node_create_binary_ordered(primary, secondary);
        if (node == NULL){
            reb_integrator_whfast_hj_node_free(primary);
            reb_integrator_whfast_hj_node_free(secondary);
            reb_integrator_whfast_hj_tree_parser_error(parser, "WHFast HJ was not able to allocate memory for the tree.");
            return NULL;
        }
        return node;
    }

    if (isdigit((unsigned char)*parser->cursor)){
        unsigned long long particle_number = 0;
        while (isdigit((unsigned char)*parser->cursor)){
            const unsigned int digit = (unsigned int)(*parser->cursor - '0');
            if (particle_number > (ULLONG_MAX - digit)/10ULL){
                reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: particle index is too large.");
                return NULL;
            }
            particle_number = 10ULL*particle_number + digit;
            parser->cursor++;
        }

        if (particle_number == 0ULL || particle_number > (unsigned long long)parser->r->N){
            reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: particle index out of range.");
            return NULL;
        }

        const size_t particle_index = (size_t)(particle_number - 1ULL);
        if (parser->used[particle_index]){
            reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: duplicate particle index.");
            return NULL;
        }

        parser->used[particle_index] = 1;
        parser->leaf_count++;

        struct hj_node* const node = reb_integrator_whfast_hj_node_create_leaf(parser->r->particles[particle_index], (int)particle_index);
        if (node == NULL){
            reb_integrator_whfast_hj_tree_parser_error(parser, "WHFast HJ was not able to allocate memory for the tree.");
            return NULL;
        }
        return node;
    }

    reb_integrator_whfast_hj_tree_parser_error(parser, "Invalid WHFast HJ tree: expected '[' or particle index.");
    return NULL;
}

REB_API int reb_integrator_whfast_hj_set_tree(struct reb_simulation* const r, const char* const tree)
{
    if (r == NULL || tree == NULL){
        return 1;
    }
    if (r->integrator.name == NULL || strcmp(r->integrator.name, "whfast_hj") != 0 || r->integrator.state == NULL){
        reb_simulation_error(r, "WHFast HJ tree can only be set when the selected integrator is whfast_hj.");
        return 1;
    }
    if (r->N > (size_t)INT_MAX + 1U){
        reb_simulation_error(r, "WHFast HJ tree does not support this many particles.");
        return 1;
    }

    unsigned char* const used = (r->N > 0) ? calloc(r->N, sizeof(unsigned char)) : NULL;
    if (r->N > 0 && used == NULL){
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree validation.");
        return 1;
    }

    struct reb_integrator_whfast_hj_tree_parser parser = {
        .r = r,
        .cursor = tree,
        .used = used,
        .leaf_count = 0,
        .error = NULL,
    };

    struct hj_node* const root = reb_integrator_whfast_hj_parse_tree_node(&parser);
    if (parser.error == NULL){
        reb_integrator_whfast_hj_tree_parser_skip_space(&parser);
        if (*parser.cursor != '\0'){
            reb_integrator_whfast_hj_tree_parser_error(&parser, "Invalid WHFast HJ tree: trailing characters.");
        }else if (root == NULL && r->N != 0){
            reb_integrator_whfast_hj_tree_parser_error(&parser, "Invalid WHFast HJ tree: empty tree for a non-empty simulation.");
        }else if (parser.leaf_count != r->N){
            reb_integrator_whfast_hj_tree_parser_error(&parser, "Invalid WHFast HJ tree: tree must include every particle exactly once.");
        }
    }

    free(used);

    if (parser.error != NULL){
        reb_integrator_whfast_hj_node_free(root);
        reb_simulation_error(r, parser.error);
        return 1;
    }

    struct reb_integrator_whfast_hj_state* const whfast = r->integrator.state;
    reb_integrator_whfast_hj_node_free(whfast->root);
    whfast->root = root;
    whfast->tree_N = r->N;
    whfast->given_tree = 1;
    return 0;
}

REB_API int reb_integrator_whfast_hj_set_binary_plus_particles_tree(struct reb_simulation* const r)
{
    if (r == NULL){
        return 1;
    }
    if (r->integrator.name == NULL || strcmp(r->integrator.name, "whfast_hj") != 0 || r->integrator.state == NULL){
        reb_simulation_error(r, "WHFast HJ tree can only be set when the selected integrator is whfast_hj.");
        return 1;
    }
    if (r->N > (size_t)INT_MAX + 1U){
        reb_simulation_error(r, "WHFast HJ tree does not support this many particles.");
        return 1;
    }

    struct hj_node* root = NULL;
    for (size_t i=0; i<r->N; i++){
        struct hj_node* const leaf = reb_integrator_whfast_hj_node_create_leaf(r->particles[i], (int)i);
        if (leaf == NULL){
            reb_integrator_whfast_hj_node_free(root);
            reb_simulation_error(r, "WHFast HJ was not able to allocate memory for the tree.");
            return 1;
        }

        if (root == NULL){
            root = leaf;
            continue;
        }

        struct hj_node* const merged = reb_integrator_whfast_hj_node_create_binary_ordered(root, leaf);
        if (merged == NULL){
            reb_integrator_whfast_hj_node_free(root);
            reb_integrator_whfast_hj_node_free(leaf);
            reb_simulation_error(r, "WHFast HJ was not able to allocate memory for the tree.");
            return 1;
        }
        root = merged;
    }

    struct reb_integrator_whfast_hj_state* const whfast = r->integrator.state;
    reb_integrator_whfast_hj_node_free(whfast->root);
    whfast->root = root;
    whfast->tree_N = r->N;
    whfast->given_tree = 1;
    return 0;
}

REB_API void reb_integrator_whfast_hj_clear_tree(struct reb_simulation* const r)
{
    if (r == NULL || r->integrator.name == NULL || strcmp(r->integrator.name, "whfast_hj") != 0 || r->integrator.state == NULL){
        return;
    }

    struct reb_integrator_whfast_hj_state* const whfast = r->integrator.state;
    reb_integrator_whfast_hj_node_free(whfast->root);
    whfast->root = NULL;
    whfast->tree_N = 0;
    whfast->given_tree = 0;
}

void reb_integrator_whfast_hj_from_inertial(struct reb_simulation* const r, struct hj_node* const node)
{
    if (node == NULL){
        return;
    }

    size_t stack_size = 0;
    size_t stack_capacity = 128;
    struct reb_integrator_whfast_hj_visit_stack_item* stack = malloc(stack_capacity*sizeof(struct reb_integrator_whfast_hj_visit_stack_item));
    if (stack == NULL){
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree traversal.");
        return;
    }
    if (reb_integrator_whfast_hj_visit_stack_push(&stack, &stack_size, &stack_capacity, node, 0)){
        free(stack);
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree traversal.");
        return;
    }

    while (stack_size > 0){
        const struct reb_integrator_whfast_hj_visit_stack_item item = stack[--stack_size];
        struct hj_node* const current = item.node;

        if (reb_integrator_whfast_hj_node_is_leaf(current)){
            if (current->particle_index >= 0 && (size_t)current->particle_index < r->N){
                current->barycenter_particle = r->particles[current->particle_index];
            }
            current->jacobi_particle = (struct reb_particle){0};
            continue;
        }

        if (!item.visited){
            if (reb_integrator_whfast_hj_visit_stack_push(&stack, &stack_size, &stack_capacity, current, 1)
                || reb_integrator_whfast_hj_visit_stack_push(&stack, &stack_size, &stack_capacity, current->secondary, 0)
                || reb_integrator_whfast_hj_visit_stack_push(&stack, &stack_size, &stack_capacity, current->primary, 0)){
                free(stack);
                reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree traversal.");
                return;
            }
            continue;
        }

        const struct reb_particle primary = current->primary->barycenter_particle;
        const struct reb_particle secondary = current->secondary->barycenter_particle;
        const double M_primary = primary.m;
        const double M_secondary = secondary.m;
        const double M_total = M_primary + M_secondary;

        if (M_total > 0.){
            current->barycenter_particle = reb_particle_com_of_pair(primary, secondary);
        }else{
            current->barycenter_particle = primary;
            current->barycenter_particle.m = 0.;
        }

        current->jacobi_particle = (struct reb_particle){0};
        current->jacobi_particle.m = (M_total > 0.) ? M_primary*M_secondary/M_total : 0.;
        current->jacobi_particle.x = secondary.x - primary.x;
        current->jacobi_particle.y = secondary.y - primary.y;
        current->jacobi_particle.z = secondary.z - primary.z;
        current->jacobi_particle.vx = secondary.vx - primary.vx;
        current->jacobi_particle.vy = secondary.vy - primary.vy;
        current->jacobi_particle.vz = secondary.vz - primary.vz;
        current->jacobi_particle.ax = secondary.ax - primary.ax;
        current->jacobi_particle.ay = secondary.ay - primary.ay;
        current->jacobi_particle.az = secondary.az - primary.az;
    }

    free(stack);
}

static void reb_integrator_whfast_hj_reconstruct_child_barycenters(struct hj_node* const node)
{
    struct reb_particle primary = node->primary->barycenter_particle;
    struct reb_particle secondary = node->secondary->barycenter_particle;
    const struct reb_particle barycenter = node->barycenter_particle;
    const struct reb_particle jacobi = node->jacobi_particle;

    const double M_primary = primary.m;
    const double M_secondary = secondary.m;
    const double M_total = M_primary + M_secondary;

    if (M_total > 0.){
        const double primary_offset = M_secondary/M_total;
        const double secondary_offset = M_primary/M_total;

        primary.x = barycenter.x - primary_offset*jacobi.x;
        primary.y = barycenter.y - primary_offset*jacobi.y;
        primary.z = barycenter.z - primary_offset*jacobi.z;
        primary.vx = barycenter.vx - primary_offset*jacobi.vx;
        primary.vy = barycenter.vy - primary_offset*jacobi.vy;
        primary.vz = barycenter.vz - primary_offset*jacobi.vz;
        primary.ax = barycenter.ax - primary_offset*jacobi.ax;
        primary.ay = barycenter.ay - primary_offset*jacobi.ay;
        primary.az = barycenter.az - primary_offset*jacobi.az;

        secondary.x = barycenter.x + secondary_offset*jacobi.x;
        secondary.y = barycenter.y + secondary_offset*jacobi.y;
        secondary.z = barycenter.z + secondary_offset*jacobi.z;
        secondary.vx = barycenter.vx + secondary_offset*jacobi.vx;
        secondary.vy = barycenter.vy + secondary_offset*jacobi.vy;
        secondary.vz = barycenter.vz + secondary_offset*jacobi.vz;
        secondary.ax = barycenter.ax + secondary_offset*jacobi.ax;
        secondary.ay = barycenter.ay + secondary_offset*jacobi.ay;
        secondary.az = barycenter.az + secondary_offset*jacobi.az;
    }else{
        primary.x = barycenter.x;
        primary.y = barycenter.y;
        primary.z = barycenter.z;
        primary.vx = barycenter.vx;
        primary.vy = barycenter.vy;
        primary.vz = barycenter.vz;
        primary.ax = barycenter.ax;
        primary.ay = barycenter.ay;
        primary.az = barycenter.az;

        secondary.x = barycenter.x + jacobi.x;
        secondary.y = barycenter.y + jacobi.y;
        secondary.z = barycenter.z + jacobi.z;
        secondary.vx = barycenter.vx + jacobi.vx;
        secondary.vy = barycenter.vy + jacobi.vy;
        secondary.vz = barycenter.vz + jacobi.vz;
        secondary.ax = barycenter.ax + jacobi.ax;
        secondary.ay = barycenter.ay + jacobi.ay;
        secondary.az = barycenter.az + jacobi.az;
    }

    node->primary->barycenter_particle = primary;
    node->secondary->barycenter_particle = secondary;
}

void reb_integrator_whfast_hj_to_inertial(struct reb_simulation* const r, struct hj_node* const node)
{
    if (node == NULL){
        return;
    }

    size_t stack_size = 0;
    size_t stack_capacity = 128;
    struct hj_node** stack = malloc(stack_capacity*sizeof(struct hj_node*));
    if (stack == NULL){
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree traversal.");
        return;
    }
    if (reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, node)){
        free(stack);
        reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree traversal.");
        return;
    }

    while (stack_size > 0){
        struct hj_node* const current = stack[--stack_size];
        if (reb_integrator_whfast_hj_node_is_leaf(current)){
            if (current->particle_index >= 0 && (size_t)current->particle_index < r->N){
                struct reb_particle* const particle = &r->particles[current->particle_index];
                const struct reb_particle source = current->barycenter_particle;
                particle->x = source.x;
                particle->y = source.y;
                particle->z = source.z;
                particle->vx = source.vx;
                particle->vy = source.vy;
                particle->vz = source.vz;
            }
            continue;
        }

        reb_integrator_whfast_hj_reconstruct_child_barycenters(current);
        if (reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, current->secondary)
            || reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, current->primary)){
            free(stack);
            reb_simulation_error(r, "WHFast HJ was not able to allocate memory for tree traversal.");
            return;
        }
    }

    free(stack);
}

/***************************** 
 * Interaction Hamiltonian  */
static void reb_integrator_whfast_hj_interaction_step_node(const struct reb_simulation* const r, struct hj_node* const node, const double _dt){
    if (node == NULL || reb_integrator_whfast_hj_node_is_leaf(node)){
        return;
    }

    size_t stack_size = 0;
    size_t stack_capacity = 128;
    struct hj_node** stack = malloc(stack_capacity*sizeof(struct hj_node*));
    if (stack == NULL){
        return;
    }
    if (reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, node)){
        free(stack);
        return;
    }

    while (stack_size > 0){
        struct hj_node* const current = stack[--stack_size];
        if (reb_integrator_whfast_hj_node_is_leaf(current)){
            continue;
        }

        struct reb_particle* const p = &current->jacobi_particle;
        const double eta = current->primary->barycenter_particle.m + current->secondary->barycenter_particle.m;
        p->vx += _dt*p->ax;
        p->vy += _dt*p->ay;
        p->vz += _dt*p->az;

        const double rj2i = 1./(p->x*p->x + p->y*p->y + p->z*p->z);
        const double rji = sqrt(rj2i);
        const double prefac = _dt*r->G*eta*rji*rj2i;
        p->vx += prefac*p->x;
        p->vy += prefac*p->y;
        p->vz += prefac*p->z;

        if (reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, current->secondary)
            || reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, current->primary)){
            free(stack);
            return;
        }
    }

    free(stack);
}

void reb_integrator_whfast_hj_interaction_step(struct reb_simulation* const r, struct hj_node* const root, const double _dt){
    reb_integrator_whfast_hj_from_inertial(r, root);
    reb_integrator_whfast_hj_interaction_step_node(r, root, _dt);
}

/***************************** 
 * DKD Scheme                */

void reb_integrator_whfast_hj_kepler_step(const struct reb_simulation* const r, struct hj_node* const node, const double _dt){
    if (node == NULL || reb_integrator_whfast_hj_node_is_leaf(node)){
        return;
    }

    size_t stack_size = 0;
    size_t stack_capacity = 128;
    struct hj_node** stack = malloc(stack_capacity*sizeof(struct hj_node*));
    if (stack == NULL){
        return;
    }
    if (reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, node)){
        free(stack);
        return;
    }

    while (stack_size > 0){
        struct hj_node* const current = stack[--stack_size];
        if (reb_integrator_whfast_hj_node_is_leaf(current)){
            continue;
        }

        const double eta = current->primary->barycenter_particle.m + current->secondary->barycenter_particle.m;
        reb_integrator_whfast_kepler_solver(&current->jacobi_particle, eta*r->G, _dt, r);

        if (reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, current->secondary)
            || reb_integrator_whfast_hj_node_stack_push(&stack, &stack_size, &stack_capacity, current->primary)){
            free(stack);
            return;
        }
    }

    free(stack);
}

void reb_integrator_whfast_hj_com_step(const struct reb_simulation* const r, struct hj_node* const root, const double _dt){
    (void)r;
    if (root == NULL){
        return;
    }
    root->barycenter_particle.x += _dt*root->barycenter_particle.vx;
    root->barycenter_particle.y += _dt*root->barycenter_particle.vy;
    root->barycenter_particle.z += _dt*root->barycenter_particle.vz;
}
//----------------------------------------------------------------------------

void reb_integrator_whfast_hj_step(struct reb_simulation* const r, void* state){
    struct reb_integrator_whfast_hj_state* whfast = state;
    const double dt = r->dt;
    if (!whfast->given_tree){
        if (reb_integrator_whfast_hj_build_tree(r, whfast)){
            reb_simulation_error(r, "WHFast HJ was not able to allocate memory for the tree.");
            return;
        }
    }else{
        if (whfast->root == NULL && r->N > 0){
            reb_simulation_error(r, "WHFast HJ given_tree is enabled, but no tree has been set.");
            return;
        }
        if (whfast->tree_N != r->N){
            reb_simulation_error(r, "WHFast HJ given_tree is enabled, but the particle count changed after the tree was set.");
            return;
        }
    }

    reb_integrator_whfast_hj_from_inertial(r, whfast->root);

    reb_integrator_whfast_hj_kepler_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_com_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast->root);

    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    reb_simulation_update_acceleration(r);
    reb_integrator_whfast_hj_interaction_step(r, whfast->root, dt);

    reb_integrator_whfast_hj_kepler_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_com_step(r, whfast->root, dt/2.);
    reb_integrator_whfast_hj_to_inertial(r, whfast->root);

    r->t += dt;
    r->dt_last_done = dt;
}
