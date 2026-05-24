/**
 * Unit test for reb_integrator_whfast_hj_build_tree().
 *
 * Five total particles, with no close encounters between non-central
 * particles. The expected tree is root -> one direct leaf child for every
 * particle.
 *
 * Expected tree:
 *
 * root
 * |-- leaf: particle 0
 * |-- leaf: particle 1
 * |-- leaf: particle 2
 * |-- leaf: particle 3
 * `-- leaf: particle 4
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

static struct reb_simulation* setup_no_close_encounter_simulation(void){
    struct reb_simulation* const r = reb_simulation_create();
    r->G = 1.0;

    add_central_particle(r);

    const double m = 5.e-8;
    add_circular_particle(r, m, 0.00);
    add_circular_particle(r, m, 0.05);
    add_circular_particle(r, m, 0.10);
    add_circular_particle(r, m, 0.15);

    return r;
}

static int assert_leaf(const struct hj_node* const node, const int particle_index){
    if (node->N_children != 0){
        return fail("expected a leaf node");
    }
    if (node->children != NULL){
        return fail("leaf node should not have children");
    }
    if (node->particle_index != particle_index){
        return fail("leaf node has wrong particle_index");
    }
    return 1;
}

static int assert_no_close_encounter_tree(const struct hj_node* const root){
    if (root->particle_index != -1){
        return fail("root should be a wrapper node");
    }
    if (root->N_children != N_PARTICLES){
        return fail("root should have one direct leaf child for every particle");
    }
    if (root->children == NULL){
        return fail("root has children count but NULL children pointer");
    }

    for (int i=0; i<N_PARTICLES; i++){
        if (!assert_leaf(&root->children[i], i)){
            return 0;
        }
    }
    return 1;
}

static int test_build_tree_no_close_encounters(void){
    struct reb_simulation* const r = setup_no_close_encounter_simulation();
    struct reb_integrator_whfast_hj_state whfast = {0};

    int ok = 1;
    if (reb_integrator_whfast_hj_build_tree(r, &whfast) != 0){
        ok = fail("reb_integrator_whfast_hj_build_tree returned an error");
        goto cleanup;
    }

    ok = assert_no_close_encounter_tree(&whfast.root);

cleanup:
    reb_integrator_whfast_hj_node_free(&whfast.root);
    reb_simulation_free(r);
    return ok;
}

int main(void){
    if (!test_build_tree_no_close_encounters()){
        return 1;
    }
    puts("reb_integrator_whfast_hj_build_tree no-close-encounter test passed.");
    return 0;
}
