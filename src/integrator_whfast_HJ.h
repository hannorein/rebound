#ifndef _INTEGRATOR_WHFAST_HJ_H
#define _INTEGRATOR_WHFAST_HJ_H

extern const struct reb_integrator reb_integrator_whfast_hj;
REB_API int reb_integrator_whfast_hj_tree_to_string(struct reb_simulation* const r, char* const buffer, const size_t buffer_size);
REB_API int reb_integrator_whfast_hj_set_tree(struct reb_simulation* const r, const char* const tree);
REB_API int reb_integrator_whfast_hj_set_binary_plus_particles_tree(struct reb_simulation* const r);
REB_API void reb_integrator_whfast_hj_clear_tree(struct reb_simulation* const r);

struct hj_node
{
    struct reb_particle barycenter_particle;
    struct reb_particle jacobi_particle;

    struct hj_node* primary;
    struct hj_node* secondary;

    int particle_index; // -1 for internal node, >= 0 for leaf
};

struct reb_integrator_whfast_hj_state {
    // Internal use
    struct hj_node* root;
    size_t tree_N;

    // Use root as a user-supplied tree instead of rebuilding the tree each step.
    int given_tree;
};

#endif
