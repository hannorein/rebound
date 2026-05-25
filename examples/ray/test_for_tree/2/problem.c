/**
 * Unit test for reb_integrator_whfast_hj_build_tree().
 *
 * Five total particles: particle 0 is the central mass, and particles 1-4
 * are on increasingly wider circular orbits with no close encounters between
 * non-central particles.
 *
 * Expected binary tree:
 *
 * root
 * |-- primary: internal {0, 1, 2, 3}
 * |   |-- primary: internal {0, 1, 2}
 * |   |   |-- primary: internal {0, 1}
 * |   |   |   |-- primary: leaf particle 0
 * |   |   |   `-- secondary: leaf particle 1
 * |   |   `-- secondary: leaf particle 2
 * |   `-- secondary: leaf particle 3
 * `-- secondary: leaf particle 4
 */

#include <math.h>
#include <stdio.h>

#include "rebound.h"
#include "../../../../src/integrator_whfast_HJ.c"

#define N_PARTICLES 5

static int fail(const char* const message){
    fprintf(stderr, "FAIL: %s\n", message);
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

static struct reb_simulation* setup_no_close_encounter_simulation(void){
    struct reb_simulation* const r = reb_simulation_create();
    r->G = 1.0;

    add_central_particle(r);

    const double m = 5.e-8;
    add_circular_particle(r, m, 1.0, 0.0);
    add_circular_particle(r, m, 2.0, 0.0);
    add_circular_particle(r, m, 3.0, 0.0);
    add_circular_particle(r, m, 4.0, 0.0);

    return r;
}

static int assert_leaf(const struct hj_node* const node, const int particle_index){
    if (node == NULL){
        return fail("expected a leaf node but found NULL");
    }
    if (node->primary != NULL || node->secondary != NULL){
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
    if (node->primary == NULL && node->secondary == NULL){
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

static int assert_tree_invariants(const struct hj_node* const node){
    if (node == NULL){
        return fail("root should not be NULL");
    }
    if (node->primary == NULL && node->secondary == NULL){
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

static int assert_no_close_encounter_tree(const struct hj_node* const root){
    if (!assert_internal(root)){
        return 0;
    }

    const struct hj_node* const node_0123 = root->primary;
    if (!assert_internal(node_0123) || !assert_leaf(root->secondary, 4)){
        return 0;
    }

    const struct hj_node* const node_012 = node_0123->primary;
    if (!assert_internal(node_012) || !assert_leaf(node_0123->secondary, 3)){
        return 0;
    }

    const struct hj_node* const node_01 = node_012->primary;
    if (!assert_internal(node_01) || !assert_leaf(node_012->secondary, 2)){
        return 0;
    }

    if (!assert_leaf(node_01->primary, 0) || !assert_leaf(node_01->secondary, 1)){
        return 0;
    }
    return 1;
}

static int assert_jacobi_link_matches_particles(const struct hj_node* const node,
        const struct reb_particle primary, const struct reb_particle secondary){
    if (!assert_internal(node)){
        return 0;
    }
    const struct reb_particle link = node->jacobi_particle;
    if (!nearly_equal(link.m, secondary.m)
            || !nearly_equal(link.x, secondary.x - primary.x)
            || !nearly_equal(link.y, secondary.y - primary.y)
            || !nearly_equal(link.z, secondary.z - primary.z)
            || !nearly_equal(link.vx, secondary.vx - primary.vx)
            || !nearly_equal(link.vy, secondary.vy - primary.vy)
            || !nearly_equal(link.vz, secondary.vz - primary.vz)){
        return fail("local Jacobi link does not match secondary minus primary");
    }
    return 1;
}

static int test_build_tree_empty_simulation(void){
    struct reb_simulation* const r = reb_simulation_create();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error for empty simulation");
        goto cleanup;
    }
    if (whfast.root != NULL){
        ok = fail("empty simulation should leave root NULL");
    }

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

static int test_build_tree_single_particle(void){
    struct reb_simulation* const r = reb_simulation_create();
    struct reb_integrator_whfast_hj_state whfast = {0};

    add_central_particle(r);

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error for one particle");
        goto cleanup;
    }

    ok = assert_leaf(whfast.root, 0)
        && assert_tree_invariants(whfast.root)
        && assert_particles_0_to_n_minus_1_once(whfast.root, 1);

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

static int test_build_tree_no_close_encounters(void){
    struct reb_simulation* const r = setup_no_close_encounter_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error");
        goto cleanup;
    }

    ok = assert_tree_invariants(whfast.root)
        && assert_each_particle_once(whfast.root)
        && assert_no_close_encounter_tree(whfast.root);

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

static int test_tree_jacobi_roundtrip(void){
    struct reb_simulation* const r = setup_no_close_encounter_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};
    struct reb_particle original[N_PARTICLES];

    for (int i=0; i<N_PARTICLES; i++){
        original[i] = r->particles[i];
    }

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error");
        goto cleanup;
    }
    reb_integrator_whfast_hj_from_inertial(r, whfast.root);

    const struct hj_node* const node_01 = whfast.root->primary->primary->primary;
    if (!assert_jacobi_link_matches_particles(node_01, original[0], original[1])){
        ok = 0;
        goto cleanup;
    }
    if (!assert_tree_invariants(whfast.root)){
        ok = 0;
        goto cleanup;
    }

    for (int i=0; i<N_PARTICLES; i++){
        r->particles[i].x += 100. + i;
        r->particles[i].y -= 50. + i;
        r->particles[i].vx -= 10. + i;
        r->particles[i].vy += 20. + i;
    }

    reb_integrator_whfast_hj_to_inertial(r, whfast.root);
    for (int i=0; i<N_PARTICLES; i++){
        if (!assert_particle_close(r->particles[i], original[i], "round-tripped particle")){
            ok = 0;
            goto cleanup;
        }
    }

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

static int test_build_tree_twice(void){
    struct reb_simulation* const r = setup_no_close_encounter_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("first reb_integrator_whfast_hj_build_tree call returned an error");
        goto cleanup;
    }
    if (!assert_tree_invariants(whfast.root) || !assert_each_particle_once(whfast.root)){
        ok = 0;
        goto cleanup;
    }

    r->particles[4].x = 8.0;
    r->particles[4].y = 0.0;
    r->particles[4].vx = 0.0;
    r->particles[4].vy = sqrt(1.0/8.0);

    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("second reb_integrator_whfast_hj_build_tree call returned an error");
        goto cleanup;
    }
    ok = assert_tree_invariants(whfast.root)
        && assert_each_particle_once(whfast.root);

cleanup:
    reb_integrator_whfast_hj_node_free(whfast.root);
    reb_simulation_free(r);
    return ok;
}

int main(void){
    if (!test_build_tree_empty_simulation()){
        return 1;
    }
    if (!test_build_tree_single_particle()){
        return 1;
    }
    if (!test_build_tree_no_close_encounters()){
        return 1;
    }
    if (!test_tree_jacobi_roundtrip()){
        return 1;
    }
    if (!test_build_tree_twice()){
        return 1;
    }
    puts("reb_integrator_whfast_hj_build_tree count/rebuild tests passed.");
    return 0;
}
