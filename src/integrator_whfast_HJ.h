#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;

struct reb_integrator_whfast_hj_node {
    size_t primary;
    size_t secondary;
    size_t barycenter;  // Particle index, or tree_N + shared barycenter index.
    struct reb_particle q;  // Relative coordinates; m is the orbit's total mass.
};

struct reb_integrator_whfast_hj_state {
    // One entry per binary orbit, in postorder. Child indices below tree_N
    // refer to particles; larger indices refer to nodes[index - tree_N].
    char* tree;
    struct reb_integrator_whfast_hj_node* nodes;
    struct reb_particle* barycenters; // Only binaries with two nonzero-mass children.
    size_t tree_N;
    size_t N_barycenters;
};

#endif
