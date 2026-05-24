/**
 * Unit test for reb_integrator_whfast_hj_build_tree().
 *
 * This intentionally includes the implementation file so this test can call
 * the file-local tree builder without changing the production API.
 *
 * Expected tree:
 *
 * root
 * |-- leaf: particle 0
 * |-- subtree: {1, 2, 3}
 * |   |-- leaf: particle 1
 * |   `-- subtree: {2, 3}
 * |       |-- leaf: particle 2
 * |       `-- leaf: particle 3
 * |-- subtree: {4, 5}
 * |   |-- leaf: particle 4
 * |   `-- leaf: particle 5
 * `-- leaf: particle 6
 */

#include <math.h>
#include <stdio.h>

#include "rebound.h"
#include "../../../../src/integrator_whfast_HJ.c"

#define N_PARTICLES 7
#define N_COMPONENTS 3
#define BIT(i) (1u << (i))

static int fail(const char* const message){
    fprintf(stderr, "FAIL: %s\n", message);
    return 0;
}

static void add_central_particle(struct reb_simulation* const r){
    struct reb_particle p = {0};
    p.m = 1.0;
    reb_simulation_add(r, p);
}

static void add_circular_particle(struct reb_simulation* const r, const double mass, const double f){
    struct reb_particle p = {0};
    p.m = mass;
    p.x = cos(f);
    p.y = sin(f);
    p.vx = -sin(f);
    p.vy = cos(f);
    reb_simulation_add(r, p);
}

static struct reb_simulation* setup_close_component_simulation(void){
    struct reb_simulation* const r = reb_simulation_create();
    r->G = 1.0;

    add_central_particle(r);

    const double m = 5.e-8;
    add_circular_particle(r, m, 0.000);  // 1: close to 2 and 3
    add_circular_particle(r, m, 0.003);  // 2: close to 1 and 3
    add_circular_particle(r, m, -0.003); // 3: close to 1 and 2
    add_circular_particle(r, m, 0.100);  // 4: close to 5
    add_circular_particle(r, m, 0.104);  // 5: close to 4
    add_circular_particle(r, m, 0.200);  // 6: isolated from other non-central particles

    return r;
}

static int collect_leaf_counts(const struct hj_node* const node, int counts[N_PARTICLES]){
    if (node->N_children == 0){
        if (node->particle_index < 0 || node->particle_index >= N_PARTICLES){
            return fail("leaf has invalid particle_index");
        }
        counts[node->particle_index]++;
        return 1;
    }

    if (node->particle_index != -1){
        return fail("internal node should have particle_index == -1");
    }
    if (node->children == NULL){
        return fail("internal node has children count but NULL children pointer");
    }

    for (int i=0; i<node->N_children; i++){
        if (!collect_leaf_counts(&node->children[i], counts)){
            return 0;
        }
    }
    return 1;
}

static int collect_mask(const struct hj_node* const node, unsigned int* const mask){
    if (node->N_children == 0){
        if (node->particle_index <= 0 || node->particle_index >= N_PARTICLES){
            return fail("non-root component contains particle 0 or invalid leaf");
        }
        *mask |= BIT(node->particle_index);
        return 1;
    }

    if (node->particle_index != -1){
        return fail("component subtree has a non-wrapper internal node");
    }

    for (int i=0; i<node->N_children; i++){
        if (!collect_mask(&node->children[i], mask)){
            return 0;
        }
    }
    return 1;
}

static int matches_expected_component(const unsigned int mask, int matched[N_COMPONENTS]){
    const unsigned int expected[N_COMPONENTS] = {
        BIT(1) | BIT(2) | BIT(3),
        BIT(4) | BIT(5),
        BIT(6),
    };

    for (int i=0; i<N_COMPONENTS; i++){
        if (!matched[i] && mask == expected[i]){
            matched[i] = 1;
            return 1;
        }
    }
    return 0;
}

static int is_leaf_particle(const struct hj_node* const node, const int particle_index){
    return node->N_children == 0 && node->particle_index == particle_index;
}

static int assert_subtree_23_shape(const struct hj_node* const node){
    int has_leaf_2 = 0;
    int has_leaf_3 = 0;

    if (node->particle_index != -1){
        return fail("subtree {2,3} should be a wrapper node");
    }
    if (node->N_children != 2){
        return fail("subtree {2,3} should have exactly two leaf children");
    }
    if (node->children == NULL){
        return fail("subtree {2,3} has NULL children pointer");
    }

    for (int i=0; i<node->N_children; i++){
        if (is_leaf_particle(&node->children[i], 2)){
            has_leaf_2 = 1;
        }else if (is_leaf_particle(&node->children[i], 3)){
            has_leaf_3 = 1;
        }else{
            return fail("subtree {2,3} should contain only leaf particles 2 and 3");
        }
    }

    if (!has_leaf_2 || !has_leaf_3){
        return fail("particles 2 and 3 should be leaf nodes inside subtree {2,3}");
    }
    return 1;
}

static int assert_123_component_shape(const struct hj_node* const node){
    int has_leaf_1 = 0;
    int has_subtree_23 = 0;

    for (int i=0; i<node->N_children; i++){
        const struct hj_node* const child = &node->children[i];
        if (is_leaf_particle(child, 1)){
            has_leaf_1 = 1;
            continue;
        }

        unsigned int mask = 0;
        if (!collect_mask(child, &mask)){
            return 0;
        }
        if (mask == (BIT(2) | BIT(3))){
            if (!assert_subtree_23_shape(child)){
                return 0;
            }
            has_subtree_23 = 1;
        }else{
            return fail("component {1,2,3} should contain leaf 1 and subtree {2,3}");
        }
    }

    if (!has_leaf_1){
        return fail("component {1,2,3} should have particle 1 as a direct leaf");
    }
    if (!has_subtree_23){
        return fail("particles 2 and 3 should form subtree {2,3}");
    }
    return 1;
}

static int assert_45_component_shape(const struct hj_node* const node){
    int has_leaf_4 = 0;
    int has_leaf_5 = 0;

    if (node->particle_index != -1){
        return fail("component {4,5} should be a wrapper node");
    }
    if (node->children == NULL){
        return fail("component {4,5} has NULL children pointer");
    }

    for (int i=0; i<node->N_children; i++){
        if (is_leaf_particle(&node->children[i], 4)){
            has_leaf_4 = 1;
        }else if (is_leaf_particle(&node->children[i], 5)){
            has_leaf_5 = 1;
        }else{
            return fail("component {4,5} should contain only leaf particles 4 and 5");
        }
    }

    if (!has_leaf_4 || !has_leaf_5){
        return fail("particle 5 should be a leaf node in component {4,5}");
    }
    return 1;
}

static int assert_root_shape(const struct hj_node* const root){
    if (root->particle_index != -1){
        return fail("root should be a wrapper node");
    }
    if (root->N_children != N_COMPONENTS + 1){
        return fail("root should have central leaf plus one child per non-central component");
    }
    if (root->children == NULL){
        return fail("root has children count but NULL children pointer");
    }

    const struct hj_node* const central = &root->children[0];
    if (central->N_children != 0 || central->particle_index != 0){
        return fail("root first child should be leaf particle 0");
    }

    int counts[N_PARTICLES] = {0};
    if (!collect_leaf_counts(root, counts)){
        return 0;
    }
    for (int i=0; i<N_PARTICLES; i++){
        if (counts[i] != 1){
            return fail("each particle should appear exactly once in the tree");
        }
    }

    int matched[N_COMPONENTS] = {0};
    const struct hj_node* component_123 = NULL;
    const struct hj_node* component_45 = NULL;
    const struct hj_node* component_6 = NULL;
    for (int i=1; i<root->N_children; i++){
        unsigned int mask = 0;
        if (!collect_mask(&root->children[i], &mask)){
            return 0;
        }
        if (mask == (BIT(1) | BIT(2) | BIT(3))){
            component_123 = &root->children[i];
        }else if (mask == (BIT(4) | BIT(5))){
            component_45 = &root->children[i];
        }else if (mask == BIT(6)){
            component_6 = &root->children[i];
        }
        if (!matches_expected_component(mask, matched)){
            return fail("root child does not match an expected close-encounter component");
        }
    }

    for (int i=0; i<N_COMPONENTS; i++){
        if (!matched[i]){
            return fail("missing expected close-encounter component");
        }
    }
    if (component_123 == NULL){
        return fail("missing component {1,2,3}");
    }
    if (component_45 == NULL){
        return fail("missing component {4,5}");
    }
    if (component_6 == NULL){
        return fail("missing component {6}");
    }
    if (!assert_123_component_shape(component_123)){
        return 0;
    }
    if (!assert_45_component_shape(component_45)){
        return 0;
    }
    if (!is_leaf_particle(component_6, 6)){
        return fail("particle 6 should be a direct leaf child of root");
    }
    return 1;
}

static int test_build_tree_components(void){
    struct reb_simulation* const r = setup_close_component_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error");
        goto cleanup;
    }

    ok = assert_root_shape(&whfast.root);

cleanup:
    reb_integrator_whfast_hj_node_free(&whfast.root);
    reb_simulation_free(r);
    return ok;
}

int main(void){
    if (!test_build_tree_components()){
        return 1;
    }
    puts("reb_integrator_whfast_hj_build_tree tree test passed.");
    return 0;
}
