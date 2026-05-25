/**
 * Unit test for reb_integrator_whfast_hj_build_tree().
 *
 * This intentionally includes the implementation file so this test can call
 * the file-local tree builder without changing the production API.
 *
 * The setup has two short-period local binaries/components:
 *
 * - particles 1, 2, and 3 form the innermost short-period component
 * - particles 4 and 5 form a wider short-period pair
 * - particle 6 is the widest isolated particle
 *
 * Expected binary tree:
 *
 * root
 * |-- primary: internal {0, 1, 2, 3, 4, 5}
 * |   |-- primary: internal {0, 1, 2, 3}
 * |   |   |-- primary: leaf particle 0
 * |   |   `-- secondary: internal {1, 2, 3}
 * |   |       |-- primary: internal {1, 2}
 * |   |       `-- secondary: leaf particle 3
 * |   `-- secondary: internal {4, 5}
 * `-- secondary: leaf particle 6
 */

#include <math.h>
#include <stdio.h>

#include "rebound.h"
#include "../../../../src/integrator_whfast_HJ.c"

#define N_PARTICLES 7
#define BIT(i) (1u << (i))
#define MASK_12 (BIT(1) | BIT(2))
#define MASK_123 (MASK_12 | BIT(3))
#define MASK_45 (BIT(4) | BIT(5))
#define MASK_0123 (BIT(0) | MASK_123)
#define MASK_012345 (MASK_0123 | MASK_45)
#define MASK_ALL (MASK_012345 | BIT(6))

static int fail(const char* const message){
    fprintf(stderr, "FAIL: %s\n", message);
    return 0;
}

static int fail_mask(const char* const label, const unsigned int actual, const unsigned int expected){
    fprintf(stderr, "FAIL: %s has mask 0x%x, expected 0x%x\n", label, actual, expected);
    return 0;
}

static int nearly_equal(const double a, const double b){
    const double scale = 1. + fabs(a) + fabs(b);
    return fabs(a - b) <= 1.e-12*scale;
}

static int assert_particle_close(const struct reb_particle actual, const struct reb_particle expected, const char* const label){
    if (!nearly_equal(actual.m, expected.m)
            || !nearly_equal(actual.x, expected.x)
            || !nearly_equal(actual.y, expected.y)
            || !nearly_equal(actual.z, expected.z)
            || !nearly_equal(actual.vx, expected.vx)
            || !nearly_equal(actual.vy, expected.vy)
            || !nearly_equal(actual.vz, expected.vz)
            || !nearly_equal(actual.ax, expected.ax)
            || !nearly_equal(actual.ay, expected.ay)
            || !nearly_equal(actual.az, expected.az)){
        fprintf(stderr, "FAIL: %s barycenter does not match its children\n", label);
        return 0;
    }
    return 1;
}

static void add_central_particle(struct reb_simulation* const r){
    struct reb_particle p = {0};
    p.m = 1.0;
    reb_simulation_add(r, p);
}

static void add_circular_particle(struct reb_simulation* const r, const double mass, const double radius, const double f){
    struct reb_particle p = {0};
    p.m = mass;
    p.x = radius*cos(f);
    p.y = radius*sin(f);
    const double speed = sqrt(1.0/radius);
    p.vx = -speed*sin(f);
    p.vy = speed*cos(f);
    reb_simulation_add(r, p);
}

static struct reb_simulation* setup_close_component_simulation(void){
    struct reb_simulation* const r = reb_simulation_create();
    r->G = 1.0;

    add_central_particle(r);

    const double m = 5.e-8;
    add_circular_particle(r, m, 1.0, 0.000);  // 1: close to 2 and 3
    add_circular_particle(r, m, 1.0, 0.003);  // 2: close to 1
    add_circular_particle(r, m, 1.0, -0.003); // 3: close to component {1,2}
    add_circular_particle(r, m, 2.0, 0.100);  // 4: close to 5
    add_circular_particle(r, m, 2.0, 0.101);  // 5: close to 4
    add_circular_particle(r, m, 4.0, 0.200);  // 6: isolated widest particle

    return r;
}

static int is_leaf_node(const struct hj_node* const node){
    return node != NULL && node->primary == NULL && node->secondary == NULL;
}

static int is_leaf_particle(const struct hj_node* const node, const int particle_index){
    return is_leaf_node(node) && node->particle_index == particle_index;
}

static int assert_leaf(const struct hj_node* const node, const int particle_index){
    if (node == NULL){
        return fail("expected a leaf node but found NULL");
    }
    if (!is_leaf_node(node)){
        return fail("expected a leaf node");
    }
    if (node->particle_index != particle_index){
        return fail("leaf node has wrong particle_index");
    }
    return 1;
}

static int assert_internal(const struct hj_node* const node){
    if (node == NULL){
        return fail("expected an internal node but found NULL");
    }
    if (node->particle_index != -1){
        return fail("internal node should have particle_index == -1");
    }
    if (node->primary == NULL || node->secondary == NULL){
        return fail("internal node should have primary and secondary children");
    }
    return 1;
}

static int collect_leaf_counts(const struct hj_node* const node, int counts[N_PARTICLES]){
    if (node == NULL){
        return fail("encountered NULL tree node");
    }
    if (is_leaf_node(node)){
        if (node->particle_index < 0 || node->particle_index >= N_PARTICLES){
            return fail("leaf has invalid particle_index");
        }
        counts[node->particle_index]++;
        return 1;
    }
    if (!assert_internal(node)){
        return 0;
    }
    return collect_leaf_counts(node->primary, counts)
        && collect_leaf_counts(node->secondary, counts);
}

static int collect_mask(const struct hj_node* const node, unsigned int* const mask){
    if (node == NULL){
        return fail("encountered NULL tree node");
    }
    if (is_leaf_node(node)){
        if (node->particle_index < 0 || node->particle_index >= N_PARTICLES){
            return fail("leaf has invalid particle_index");
        }
        *mask |= BIT(node->particle_index);
        return 1;
    }
    if (!assert_internal(node)){
        return 0;
    }
    return collect_mask(node->primary, mask)
        && collect_mask(node->secondary, mask);
}

static int assert_mask(const struct hj_node* const node, const unsigned int expected, const char* const label){
    unsigned int actual = 0;
    if (!collect_mask(node, &actual)){
        return 0;
    }
    if (actual != expected){
        return fail_mask(label, actual, expected);
    }
    return 1;
}

static int assert_leaf_pair(const struct hj_node* const node, const int a, const int b, const char* const label){
    if (!assert_internal(node)){
        return 0;
    }
    if ((is_leaf_particle(node->primary, a) && is_leaf_particle(node->secondary, b))
            || (is_leaf_particle(node->primary, b) && is_leaf_particle(node->secondary, a))){
        return 1;
    }
    fprintf(stderr, "FAIL: %s should contain only leaf particles %d and %d\n", label, a, b);
    return 0;
}

static int assert_tree_invariants(const struct hj_node* const node){
    if (node == NULL){
        return fail("root should not be NULL");
    }
    if (is_leaf_node(node)){
        if (node->particle_index < 0 || node->particle_index >= N_PARTICLES){
            return fail("leaf has invalid particle_index");
        }
        return 1;
    }
    if (!assert_internal(node)){
        return 0;
    }
    if (node->primary->barycenter_particle.m + 1.e-15 < node->secondary->barycenter_particle.m){
        return fail("primary child should be at least as massive as secondary child");
    }

    const struct reb_particle expected = reb_particle_com_of_pair(
        node->primary->barycenter_particle,
        node->secondary->barycenter_particle
    );
    if (!assert_particle_close(node->barycenter_particle, expected, "internal node")){
        return 0;
    }
    return assert_tree_invariants(node->primary)
        && assert_tree_invariants(node->secondary);
}

static int assert_each_particle_once(const struct hj_node* const root){
    int counts[N_PARTICLES] = {0};
    if (!collect_leaf_counts(root, counts)){
        return 0;
    }
    for (int i=0; i<N_PARTICLES; i++){
        if (counts[i] != 1){
            return fail("each particle should appear exactly once in the tree");
        }
    }
    return 1;
}

static int assert_particles_0_to_n_minus_1_once(const struct hj_node* const root, const int N_expected){
    int counts[N_PARTICLES] = {0};
    if (N_expected < 0 || N_expected > N_PARTICLES){
        return fail("invalid expected particle count");
    }
    if (!collect_leaf_counts(root, counts)){
        return 0;
    }
    for (int i=0; i<N_PARTICLES; i++){
        const int expected_count = i < N_expected ? 1 : 0;
        if (counts[i] != expected_count){
            return fail("tree does not contain the expected particle index set");
        }
    }
    return 1;
}

static int assert_close_component_tree(const struct hj_node* const root){
    if (!assert_internal(root) || !assert_mask(root, MASK_ALL, "root")){
        return 0;
    }
    if (!assert_mask(root->primary, MASK_012345, "root primary")
            || !assert_leaf(root->secondary, 6)){
        return 0;
    }

    const struct hj_node* const node_012345 = root->primary;
    if (!assert_mask(node_012345->primary, MASK_0123, "subtree {0,1,2,3}")
            || !assert_mask(node_012345->secondary, MASK_45, "subtree {4,5}")){
        return 0;
    }

    const struct hj_node* const node_0123 = node_012345->primary;
    if (!assert_leaf(node_0123->primary, 0)
            || !assert_mask(node_0123->secondary, MASK_123, "subtree {1,2,3}")){
        return 0;
    }

    const struct hj_node* const node_123 = node_0123->secondary;
    if (!assert_mask(node_123->primary, MASK_12, "subtree {1,2}")
            || !assert_leaf(node_123->secondary, 3)){
        return 0;
    }

    return assert_leaf_pair(node_123->primary, 1, 2, "subtree {1,2}")
        && assert_leaf_pair(node_012345->secondary, 4, 5, "subtree {4,5}");
}

static struct reb_simulation* setup_unequal_mass_pair_simulation(void){
    struct reb_simulation* const r = reb_simulation_create();
    r->G = 1.0;

    add_central_particle(r);
    add_circular_particle(r, 1.e-8, 1.0, 0.000);
    add_circular_particle(r, 4.e-8, 1.0, 0.001);

    return r;
}

static int test_build_tree_components(void){
    struct reb_simulation* const r = setup_close_component_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error");
        goto cleanup;
    }

    ok = assert_tree_invariants(whfast.root)
        && assert_each_particle_once(whfast.root)
        && assert_close_component_tree(whfast.root);

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

static int test_build_tree_unequal_mass_ordering(void){
    struct reb_simulation* const r = setup_unequal_mass_pair_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error");
        goto cleanup;
    }

    ok = assert_tree_invariants(whfast.root)
        && assert_particles_0_to_n_minus_1_once(whfast.root, 3)
        && assert_internal(whfast.root)
        && assert_leaf(whfast.root->primary, 0)
        && assert_mask(whfast.root->secondary, MASK_12, "unequal-mass pair")
        && assert_leaf(whfast.root->secondary->primary, 2)
        && assert_leaf(whfast.root->secondary->secondary, 1);

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

int main(void){
    if (!test_build_tree_components()){
        return 1;
    }
    if (!test_build_tree_unequal_mass_ordering()){
        return 1;
    }
    puts("reb_integrator_whfast_hj_build_tree close-component binary tests passed.");
    return 0;
}
